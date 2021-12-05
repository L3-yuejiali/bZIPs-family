**Final Project on the query sequence of LOC5519226**
1. **use a BLAST search to find homologs of the query protein.**
uncompress the proteomes and concatenate proteomes into a single file:
```
gunzip proteomes/*.gz
cat  proteomes/* > allprotein.fas
```
build a blast database
`makeblastdb -in allprotein.fas -parse_seqids -dbtype prot`

download the query sequence form ncbi
`ncbi-acc-download -F fasta -m protein XP_032219950.1`

perform a blast search by the following command and save output into the file `Nvectensis.blastp.typical.out
`

` blastp -db allprotein.fas -query XP_032219950.1.fa -outfmt 0 -max_hsps 1 -out Nvectensis.blastp.typical.out`
 
 a detailed analysis of high scoring pairs (top scoting hit) and E-values can be processed in another output
 
` blastp -db allprotein.fas -query XP_032219950.1.fa -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle"  -max_hsps 1 -out vectensis.blastp.detail.out`

**Filter the BLAST output by setting e-value as 1e-10**
`awk '{if ($6<0.0000000001)print $1 }' Nvectensis.blastp.detail.out > Nvectensis.blastp.detail.filtered.out`
`wc -l Nvectensis.blastp.detail.filtered.out`
`seqkit grep --pattern-file Nvectensis.blastp.detail.filtered.out allprotein.fas > Nvectensis.blastp.detail.filtered.fas`

**generate multiple sequence alignment in muscle:**
`muscle -in Nvectensis.blastp.detail.filtered.fas -out Nvectensis.blastp.detail.filtered.aligned.fas`
 
statistical analysis of average percent identity by t_coffee:
`t_coffee -other_pg seq_reformat -in Nvectensis.blastp.detail.filtered.aligned.fas -output sim`
 alignment output in alv and review:
` alv -kli --majority Nvectensis.blastp.detail.filtered.aligned.fas | less -RS`
**Removed highly gapped positions in the alignment by t_coffee:**
`t_coffee -other_pg seq_reformat -in Nvectensis.blastp.detail.filtered.aligned.fas -action +rm_gap 50 -out Nvectensis.blastp.detail.filtered.aligned.r50.fas`

output:
`alv -kli --majority Nvectensis.blastp.detail.filtered.aligned.r50.fas | less -RS`

**2. Phylogenetic tree made**
Download a software **newick utilities** for visualizaing and manipulating phylogeneric trees.
```
 cd ~/tools
 git clone git://github.com/tjunier/newick_utils.git
    cd newick_utils/
     autoreconf -fi
        ./configure
            make
                sudo make install
```
test the installed programs:
`  nw_display -h`

**Constrcut a phylogenetic tree from sequences by IQ-TREE:**
`sed "s/ /_/g" Nvectensis.blastp.detail.filtered.aligned.fas > Nvectensis.blastp.detail.filtered.aligned_.fas`
`iqtree -s Nvectensis.blastp.detail.filtered.aligned_.fas -nt 2`

output in:
`Nvectensis.blastp.detail.filtered.aligned_.fas.iqtree`

**Root the optimal phylogeny using midpoint rooting method:**

`gotree reroot midpoint -i Nvectensis.blastp.detail.filtered.aligned_.fas.treefile -o Nvectensis.blastp.detail.filtered.aligned_.fas.midpoint.treefile`

output display:
`nw_display -s  Nvectensis.blastp.detail.filtered.aligned_.fas.midpoint.treefile -w 1000 -b 'opacity:0' >  Nvectensis.blastp.detail.filtered.aligned_.fas.midpoint.treefile.svg`
 
**3. reconciling a gene and species tree for evaluating the duplication and loss events:**

install some software as needed(Numpy, ETE toolkit, cargo) and thirdkind, a software for generating a graphic of tree reconciliations:
```

    sudo yum install numpy
    sudo easy_install -U ete3 
curl https://sh.rustup.rs -sSf | sh
bash
cargo install thirdkind
```
install Notung for construct genes and species reconciled tree:
`java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar --help`

reconcile the gene and species tree using Notung for 
`java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -s speciesTreeBilateriaCnidaria.tre -g Nvectensis.blastp.detail.filtered.aligned_.fas.midpoint.treefile --reconcile --speciestag prefix --savepng --events`
Generate a RecPhyloXML file to view the gene-within-species tree use the program recphylovisu:
`python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g Nvectensis.blastp.detail.filtered.aligned_.fas.midpoint.treefile.reconciled --include.species `
or use thirdkind to generate a .svg file uploaded to the repository:
`thirdkind -f Nvectensis.blastp.detail.filtered.aligned_.fas.midpoint.treefile.genes.tre.reconciled.xml -o  Nvectensis.blastp.detail.filtered.aligned_.fas.midpoint.treefile.genes.tre.reconciled.svg`
Minimizing the duplication and loss while rooting via Notung:
`java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -s speciesTreeBilateriaCnidaria.tre -g Nvectensis.blastp.detail.filtered.aligned_.fas.midpoint.treefile --root --speciestag prefix --savepng --events`

IQ-TREE re-runed to generate ultrafast bootstrap support value for optimal gene tree and re-root the tree:
`iqtree -s Nvectensis.blastp.detail.filtered.aligned_.fas -bb 1000 -nt 2 -m rtREV+F+R5 -t Nvectensis.blastp.detail.filtered.aligned_.fas.midpoint.treefile -pre Nvectensis.genes.ufboot`

`gotree reroot midpoint -i Nvectensis.genes.ufboot -o Nvectensis.genes.midpoint.ufboot`

**4. Domain prediction via Interproscan and Pfam-A** 
install DataMash utility:
```
cd ~/tools
wget http://ftp.gnu.org/gnu/datamash/datamash-1.3.tar.gz
tar -xzf datamash-1.3.tar.gz
cd datamash-1.3
./configure
make
make check
sudo make install 
```
Run iprscan5 for domain identify on each gene:
`iprscan5   --email yuejia.li@stonybrook.edu  --multifasta --useSeqId --sequence   Nvectensis.blastp.detail.filtered.fas`
concatenate all output files into a single file:
`cat ~/myfamily/*.tsv.tsv > ~/myfamily/Nvectensis.domains.all.tsv`

filter the domains only defined by Pfam database:
```
grep Pfam ~/myfamily/Nvectensis.domains.all.tsv >  ~/myfamily/Nvectensis.domains.pfam.tsv
``` 
**Plot the phylogeny with domains annotated using a web-service Evolview:**

re-arranging the interproscan output and save to the file `Nvectensis.domains.pfam.evol.tsv`(an appropriate name should be bZIPs):
`awk 'BEGIN{FS="\t"} {print $1"\t"$3"\t"$7"@"$8"@"$5}' ~/myfamily/Nvectensis.domains.pfam.tsv | datamash -sW --group=1,2 collapse 3 | sed 's/,/\t/g' | sed 's/@/,/g' > ~/myfamily/Nvectensis.domains.pfam.evol.tsv`

A rooted gene tree has been chosen from the file `maguk.blastp.detail.filtered.aligned_.fas.MendozaRoot.treefile` and annotated by uploading the file `Nvectensis.domains.pfam.evol.tsv` in Evolview

The domain-on-tree graphic plotted in Evolview and saved to the file `Nvectensis.png`

