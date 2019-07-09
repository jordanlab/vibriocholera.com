#!/home/blast/anaconda3/bin/python3
import os
import subprocess
import tempfile
from pyfaidx import Fasta
import click
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord
from bokeh.resources import CDN
from bokeh.embed import file_html
from Bio import SeqIO

@click.command()
@click.option('-i', '--input', 'input', default=None,
			  help="Input hmmer TBLOUT file", type=click.File(mode='r'))
@click.option('-f', '--fasta', 'fasta', default=None,
			  help="T6SS cluster alignment being processed, chose from [L, 1, 2, 3, 4, 5]",
			  type=str)
@click.option('-a', '--faa', 'faa', default=None,
			  help="T6SS cluster alignment being processed, chose from [L, 1, 2, 3, 4, 5]",
			  type=click.File(mode='r'))
@click.option('-n', '--fna', 'fna', default=None,
			  help="T6SS cluster alignment being processed, chose from [L, 1, 2, 3, 4, 5]",
			  type=click.File(mode='r'))
@click.option('-g', '--gff', 'gff', default=None,
			  help="T6SS cluster alignment being processed, chose from [L, 1, 2, 3, 4, 5]",
			  type=click.File(mode='r'))
@click.option('-o', '--output', 'output', default='-',
			  help="Output file path",
			  type=click.Path(exists=True))
@click.option('--name', 'name', default=None,
			  help="Output file name",
			  type=str)
@click.option('--log', 'log', default="conservation.log", type=click.File(mode="w"))
def main(input, faa, fna, fasta, gff, output, log, name):
	gff_path = gff.name
	gff_data = {}
	header = []
	mol_length = {}
	genome = Fasta(fasta)
	for rec in gff:
		rec = rec.rstrip()
		if rec.startswith("##"): continue
		elif rec.startswith("# Seq"):
			header.append(rec.rstrip())
			length = int(rec.split(";")[1].split("=")[1])
		elif rec.startswith("#"): continue
		elif rec == '"': continue
		else:
			line = rec.rstrip().split("\t")
			if line[0] not in gff_data:
				gff_data[line[0]] = {}
			if line[0] not in mol_length:
				mol_length[line[0]] = length
			id = line[8].split(";")[0].split("_")[1]
			id = f"{line[0]}_{id}"
			gff_data[line[0]][id] = line
	# print(f"Sequence length: {mol_length}")
	# print(gff_data)
	hits = {}
	prots = {}
	for line in input:
		if not line.startswith("#"):
			line = line.rstrip().split()
			prots[line[2]] = line[0]
			if line[0] in hits and line[2] not in hits[line[0]]:
				hits[line[0]].append(line[2])
			else:
				hits[line[0]] = [line[2]]
	# import pprint
	# pp = pprint.PrettyPrinter(indent=4)
	# pp.pprint(hits)
	# pp.pprint(prots)
	if "vgrg" in hits:
		for vgrg in hits['vgrg']:
			contig, id = vgrg.rsplit("_",1)
			if f"{contig}_{int(id)-1}" in prots:
				genes = build_cluster(fwd=True, start=f"{contig}_{int(id)-1}", hcp=True, annotations=hits, peptides=prots)
				# print(genes)
			elif  f"{contig}_{int(id)+1}" in prots:
				genes = build_cluster(fwd=False, start=f"{contig}_{int(id)+1}", hcp=True, annotations=hits, peptides=prots)
				# print(genes)
			else:
				genes= build_cluster(fwd=False, start=vgrg, hcp=False, annotations=hits, peptides=prots)
				# print(genes)
			if any(genes[x] in ['hydrolase', 'lipase'] for x in genes):
				print(f"Found Aux1 or Aux4")
				gff_data = add_names(genes, gff_data, "Aux 1 or 4")
			elif any(genes[x] in ['ntpase', 'transferase', 'lysm'] for x in genes):
				print(f"Found Aux2")
				gff_data = add_names(genes, gff_data, "Aux2")
			else:
				print(f"Unknown partial cluster\n{genes}")
	print("Searching for Aux3 cluster")
	cluster3 = {}
	if "aux3" in hits:
		print("Aux3 present")
		for aux3 in hits["aux3"]:
			contig, id = aux3.rsplit("_",1)
			if f"{contig}_{int(id)-1}" in prots:
				cluster3[aux3] = "effector"
				cluster3[f"{contig}_{int(id)-1}"] = "PAAR"
				cluster3[f"{contig}_{int(id)+1}"] = "immunity"
			elif f"{contig}_{int(id)+1}" in prots:
				cluster3[aux3] = "effector"
				cluster3[f"{contig}_{int(id)+1}"] = "PAAR"
				cluster3[f"{contig}_{int(id)-1}"] = "immunity"
			else:
				pass
		# print(cluster3)
		gff_data = add_names(cluster3, gff_data, "Aux3")
	else:
		print("No Aux3")

	nucl = tempfile.NamedTemporaryFile(mode="r+", suffix=".nucl")
	hmmout = tempfile.NamedTemporaryFile(mode="r+", suffix=".bed")
	aux5gff = tempfile.NamedTemporaryFile(mode="r+", suffix=".gff")
	aux5_nhmmer = f" nhmmer --cpu 4 --F1 .1 -o /dev/null --tblout {nucl.name} /home/blast/prediction_server/server/hmm_profiles//aux5 {fasta}"
	print(f"Searching for Aux5 clusters: {aux5_nhmmer}")
	subprocess.call(aux5_nhmmer, shell=True)
	aux5_pres = False
	for l in nucl:
		if l.startswith("#"): continue
		l = l.rstrip().split()
		percID = (max(int(l[7]),int(l[6])) - min(int(l[7]),int(l[6])))/6000
		if percID > 0.9:
			hmmout.write(f"{l[0]}\t{l[8]}\t{l[9]}\n")
			aux5_pres = True
	hmmout.seek(0)
	print(hmmout.read())
	bedtool_cmd = f"bedtools intersect -a {gff_path} -b {hmmout.name} > {aux5gff.name}"
	subprocess.call(bedtool_cmd, shell=True)
	if aux5_pres:
		nucleotides = Fasta(fna.name)
		aux5_fh = tempfile.NamedTemporaryFile(mode="r+", suffix=".fasta")
		for rec in aux5gff:
			rec.rstrip()
			if rec == "": continue
			else:
				line = rec.rstrip().split("\t")
				id = line[8].split(";")[0].split("_")[1]
				id = f"{line[0]}_{id}"
				aux5_fh.write(f">{id}\n{str(nucleotides[id])}\n")
		aux5_fh.seek(0)
		# print(aux5_fh.read())
		aux5hmms = ['aux5_hcp.hmm', 'aux5_vgrg.hmm', 'aux5_tap.hmm', 'aux5_eff.hmm']
		aux5hmmstat = {'aux5_hcp.hmm' : 228, 'aux5_vgrg.hmm': 2100, 'aux5_tap.hmm' : 747, 'aux5_eff.hmm': 1725}
		aux5 = {}
		for hmm in aux5hmms:
			a5hmmout = tempfile.NamedTemporaryFile(mode="r+", suffix=".out")
			a5nhmmer_cmd = f"nhmmer --F1 .1 -o /dev/null --tblout {a5hmmout.name} /home/blast/prediction_server/server/hmm_profiles/{hmm} {aux5_fh.name}"
			print(f"Running nhmmer: {a5nhmmer_cmd}")
			subprocess.call(a5nhmmer_cmd, shell=True, stderr=subprocess.DEVNULL )
			for l in a5hmmout:
				if l.startswith("#"): continue
				l = l.rstrip().split()
				percID = (max(int(l[7]),int(l[6])) - min(int(l[7]),int(l[6])))/aux5hmmstat[hmm]
				if percID > 0.9:
					print(f"Found result for {hmm}: {l[0]}")
					aux5[l[0]] = l[2].split("_")[1]
		# pp.pprint(aux5)
		gff_data = add_names(aux5, gff_data, "Aux5")
	annotation = {}
	uhOhSpaghettiOs = True
	for contig in gff_data:
		annotation[contig] = {}
		features = []
		for seq in gff_data[contig]:
			if gff_data[contig][seq][8].startswith("Name="):
				label = gff_data[contig][seq][8].split(";")[0].split("=")[1]
				# print(f"label={label},start={gff_data[contig][seq][3]}, end={gff_data[contig][seq][4]}, strand={gff_data[contig][seq][6]}1")
				features.append(
					GraphicFeature(start=int(gff_data[contig][seq][3]),
										   end=int(gff_data[contig][seq][4]),
										   strand=int(f"{gff_data[contig][seq][6]}1"),
										   color="#ffd700",
										   label=label)
					)
				annotation[contig][seq] = label
				annotation[seq] = label
			else:
				# print(f"label={None},start={gff_data[contig][seq][3]}, end={gff_data[contig][seq][4]}, strand={gff_data[contig][seq][6]}1")
					features.append(
							GraphicFeature(start=int(gff_data[contig][seq][3]),
								           end=int(gff_data[contig][seq][4]),
								           strand=int(f"{gff_data[contig][seq][6]}1")
								           )
							)

		# print(mol_length[contig])
		# sequence = str(genome[contig])
		if annotation[contig] != {}:
			uhOhSpaghettiOs = False
			record = GraphicRecord(sequence_length=mol_length[contig], features=features)
			plot = record.plot_with_bokeh(figure_width=15)
			with open(f"{output}/{contig}.html", "w+") as f:
				f.write(file_html(plot, CDN, f"{contig}"))
		proteins = Fasta(faa.name)

		with open(f"{output}/proteins.faa", "w+") as p_fh:
			for seq in proteins.keys():
				if seq in annotation:
					p_fh.write(f">{annotation[seq]}\n{str(proteins[seq])}\n")
				else:
					p_fh.write(f">{seq}\n{str(proteins[seq])}\n")
		cds = Fasta(faa.name)
		with open(f"{output}/nucleotides.fna", "w+") as n_fh:
			for seq in cds.keys():
				if seq in annotation:
					n_fh.write(f">{annotation[seq]}\n{str(cds[seq])}\n")
				else:
					n_fh.write(f">{seq}\n{str(cds[seq])}\n")


	if uhOhSpaghettiOs:
		with open(f"{output}/nohits.html", "w+") as sadface:
			sadface.write('<div class="splash"><div class="middle"><h1>No T6SS predicted :(</h1></div><div class="bottomleft"><p>T6SS.Vibriocholera.com</p></div></div>\n')


def build_cluster(fwd, start, hcp, annotations, peptides):
	contig, id = start.rsplit("_",1)
	id = int(id)
	genes = {}
	if not hcp:
		print("No HCP")
		fwd_list, rev_list = {start: peptides[start]}, {start: peptides[start]}
		count = 0
		while count <= 1:
			id += 1
			if id == 0: break
			orf = f"{contig}_{id}"
			if orf not in peptides:
				count += 1
				fwd_list[orf] = "predicted T6SS protein"
			else:
				fwd_list[orf] = peptides[orf]
			# 	count += -1
		count = 0
		while count <= 1:
			id += -1
			if id == 0: break
			orf = f"{contig}_{id}"
			if orf not in peptides:
				count += 1
				rev_list[orf] = "predicted T6SS protein"
			else:
				rev_list[orf] = peptides[orf]
			# 	count += -1
			# fwd_list.append(f"{contig}_{int(id)}")
		if len(fwd_list) > len(rev_list):
			return fwd_list
		else:
			return rev_list
	else:
		print("HCP present")
		genes[start] = peptides[start]
		count = 0
		if fwd: incr = 1
		else: incr = -1
		while count <= 1:
			id += incr
			if id == 0: break
			orf = f"{contig}_{id}"
			if orf not in peptides:
				count += 1
				genes[orf] = "predicted T6SS protein"
			else:
				genes[orf] = peptides[orf]
			# 	count += -1
			# genes.append(f"{contig}_{int(id)}")
		return genes

def add_names(genes, gff_data, cluster_name):
	for gene in genes:
		contig, id = gene.rsplit("_", 1)
		if gene in gff_data[contig]:
			info = gff_data[contig][gene][8]
			info = f"Name={gene}|{cluster_name}|{genes[gene]};{info}"
			print(info)
			gff_data[contig][gene][8] = info
	return gff_data



if __name__ == "__main__":
	main()