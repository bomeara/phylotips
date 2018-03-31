---
title: "PAUP batch files"
author: "Brian O'Meara"
date: '2018-03-30'
comment: no
contentCopyright: no
categories: paup
description: ''
flowchartDiagrams:
  enable: no
  options: ''
hiddenFromHomePage: no
keywords: []
lastmod: '2018-03-30T22:55:35-04:00'
mathjax: no
mathjaxEnableSingleDollar: no
postMetaInFooter: no
reward: no
sequenceDiagrams:
  enable: no
  options: ''
slug: paup-batch-files
tags: paup
autoCollapseToc: no
toc: no
---

This is a [backup](https://web.archive.org/web/20031220172052/http://www.brianomeara.info:80/phylogenetics.html) courtesy of the Internet Archive's [Wayback Machine](https://web.archive.org) of a page I made shortly after graduate school with batch files for PAUP. May they be of some utility to you. 

## Phylogenetic analysis

Each item includes PAUP [batch files](#batch) or commands where
appropriate. I wrote this page for the Farrell Lab website, but have
included it here as well. This page is mostly useful for the batch files
which are on it. Some useful references which provide broader overviews
of phylogenetics are *[Phylogenetic Trees Made
Easy](https://web.archive.org/web/20031220172052/http://www.sinauer.com/Titles/Text/3115.html)*
by Barry Hall, chapters 11 and 12 of the somewhat dated (1996)
*[Molecular Systematics, 2nd
Edition](https://web.archive.org/web/20031220172052/http://www.sinauer.com/Titles/frhillis.htm)*,
edited by David Hillis, Craig Moritz, and Barbara Mable, and the
recently updated [Molecular
Systematics](https://web.archive.org/web/20031220172052/http://www.bioinf.org/molsys/index.html)
website. PAUP also has a [quick start
tutorial](https://web.archive.org/web/20031220172052/http://paup.csit.fsu.edu/downl.html#Anchor-58521)
on its [web
page](https://web.archive.org/web/20031220172052/http://paup.csit.fsu.edu/)
as well as a useful [discussion
group](https://web.archive.org/web/20031220172052/http://pauptech.si.edu/%7Epaupforum/)
for asking questions. Useful software links can be found on my [general
links](/web/20031220172052/http://www.brianomeara.info:80/links.html)
page. Please direct any comments about this page to
[me](https://web.archive.org/web/20031220172052/mailto:bcomeara@ucdavis.edu).
This page was last updated on August 20, 2002 , though most of its
information has not been updated since June 10, 2002.

Each step listed individually\
[NEXUS file preparation](#nexus)\
[Batch files](#batch)\
[Testing for homogeneity, what to do if this fails](#ild)\
[Saturation plots and weighting](#saturation)\
[Exhaustive and branch-and-bound searches](#exhaustive) (\<13 taxa)\
[Heuristic searches](#heuristic)\
[Parsimony ratchet](#ratchet)\
[Bootstrap searches](#bootstrap)\
[Decay indices, including partitioned Bremer support](#decay)\
[](#bayes)[Choosing a likelihood model](phylogenetics.html#likelihood)\
[Testing for a molecular clock\
](#clock)[Calibrating using a published calibration](#Calibrate1)\
[Calibrating using fossil data](#Calibrate2)\
[Computing error bars on node times](#AgeErrors)\
[Bayesian searches](#bayes) (MrBayes)\
[Hypothesis testing](#hypothesis)\
[Tree drawing and editing](#Draw)\

[]{#nexus}NEXUS file preparation

1.  Open the NEXUS file with
    [MacClade](https://web.archive.org/web/20031220172052/http://phylogeny.arizona.edu/macclade/macclade.html)
    4.03 or higher.
2.  Select characters for one locus by shift-clicking on characters for
    that locus.
3.  Under the Characters menu, choose \"Character Set\...\" then \"Store
    Character Set\".
4.  Name this set by the name of the locus (omit spaces)
5.  Do the same for other loci.
6.  Make sure that the genetic code is correct for your loci by going to
    \"Genetic Code\...\" in the Characters menu. For mixed nuclear and
    mitochondrial loci, choose *Drosophila* mtDNA.
7.  For protein-coding loci, select all the characters in one exon, then
    go to Characters -\> \"Codon Positions\" -\> \"Calculate
    Positions\...\". Select the option \"Choose to minimize stop
    codons\".
8.  Select non-protein-coding regions (i.e., introns, rDNA) and go to
    Characters -\> \"Codon Positions\" -\> \"Non-coding\".
9.  Go to Taxa -\> \"Taxon list\".
10. Move the outgroup taxa so that they are next to each other by
    dragging the number to the left of the taxon name.
11. Select the whole group of outgroup taxa by shift-clicking on their
    numbers.
12. Go to Taxa -\> \"Taxon sets\" -\> \"Store taxon set\".
13. Name this taxon set \"outgroup\".
14. Do the same for other taxa you would like in a set. Feel free to
    change the order of the taxa to do this.
15. When you are done, click on the icon that looks like a Mayan
    pyramid, then click on one of the taxa. This alphabetizes the taxa,
    making later steps in PAUP easier.
16. Go to File -\> \"Print list\" to print the list, along with taxon
    numbers.
17. Save the file and quit MacClade.
18. Open PAUP 4.0b8 or higher. PAUP often has updates, so check its
    [website](https://web.archive.org/web/20031220172052/http://paup.csit.fsu.edu/)
    to make sure you are using a recent version. Also check for [bug
    reports](https://web.archive.org/web/20031220172052/http://paup.csit.fsu.edu/problems.html).
19. From within PAUP, open the file you just saved with MacClade. Choose
    \"Edit\", not \"Execute\".
20. Scroll down to below the data matrix.
21. NEXUS files are arranged in blocks, starting with \"BEGIN \#\#\#;\"
    and ending with \"END;\". Go to the SETS block [skip this and the
    next 2 instructions if you have only one locus or no coding
    regions].
22. Start a new line. Type \"CHARSET\", a space, then \"\#\#\#first =
    \", replacing \#\#\# with the name of your first coding locus.
23. Look for the CODONS block. You should see something within this
    block that looks like \"1: 1-1300\\3 1400-2700\\3\", with the
    numbers changed. If your first coding locus is 1302 long, for
    example, copy the \"1-1300\\3\"and place it after the equals sign in
    the line you just created. If your first coding locus extends all
    the way to 2700 in this example, just with an intron in the middle,
    copy the entire \"1-1300\\3 1400-2700\\3\" to after the equals sign
    in the line you just created.
24. Type a semicolon at the end of this line after pasting in the
    numbers.
25. You have just created a character set containing all the bases with
    the first codon position from your first coding locus. Repeat the
    same process for other loci, all three positions. The general format
    is \"CHARSET charset_name = character_numbers;\".
26. You can do the same thing with sets of taxa. You should already have
    a taxset named \"outgroup\". You can make other taxsets using the
    format \"TAXSET taxset_name = taxon_numbers;\". The taxon
    number is the same as the number on the taxon\'s left in MacClade\'s
    taxon list window, which you have printed out. It\'s also the order
    of taxa within the DATA block.
27. Create a PAUP block. Make sure this is after one of the \"end;\"
    lines but before the next \"Begin \#\#\#;\" lines. The format is
    \"Begin PAUP;\" skip a few lines, then \"end;\". Don\'t forget the
    semicolons.
28. []{#charpartition}Within the PAUP block, create a character
    partition. The general format is \"CHARPARTITION [partition name]
    = [charset name 1]:[charset name 1], [charset name
    2]:[charset name 2], [charset name 3]:[charset name 3];\" For
    example, if the data consists of 28S and COI sequence, you would
    have a charset called 28s, another charset called coI, and the
    charparition command would be \"CHARPARTITION locus = 28s:28s,
    coI:coI;\" This is used in ILD tests.
29. Within the PAUP block, paste the following text:\
    \
    outgroup outgroup /only; set maxtrees=5000 increase=auto
    autoclose=no torder=right tcompress=yes taxlabels=full
    storebrlens=yes storetreewts=yes outroot=paraphyl;\
    \
    This sets the outgroup to be whatever you determined in the outgroup
    taxset, sets a reasonable starting maxtrees number, and sets other
    PAUP preferences to useful ones for starting a search. Other PAUP
    commands can be included in the PAUP block as well.
30. You may find it useful to keep a running log of all analyses done.
    This can be helpful in checking what p-values were of various tests,
    how many replicates were done in a particular search, etc. To do
    this, paste the following text within the PAUP block:\
    \
    log file=[name] append=yes replace=no;\
    \
    If you do this, delete the log commands from other batch files (see
    below) so that everything is logged to one file.\

[]{#batch}Batch files

Batch files are ways to give PAUP commands without going through the
menus (Mac) or entering items one at a time in the command line (all
platforms). For example, you could set up a PAUP batch file which would
do a heuristic search saving the best trees from each addition sequence
replicate, save these trees to a file, filter for the best trees, and
finally save these best trees to a second tree file. Batch files are
also a good way to ensure consistency between searches \-- for example.
you won\'t have to remember if you used 5 or 10 addition sequence
replicates per bootstrap replicate, as that information will be in the
batch file. A batch file has the following format:

\#nexus [Tells PAUP that this is a NEXUS file]\
begin PAUP; [Starts the PAUP block]\
*[PAUP commands]*; [Gives PAUP instructions on what to do \-- the
same thing you enter in the command line]\
end; [Ends the PAUP block]

To create a batch file, first create a new document within PAUP or in
some other text editor (such as SimpleText). Then type \"\#nexus\" in
the first line and \"begin paup;\" in the second. After that, include
the PAUP commands you want \-- PAUP\'s command reference document,
available from the
[downloads](https://web.archive.org/web/20031220172052/http://paup.csit.fsu.edu/downl.html)
section of its
[website](https://web.archive.org/web/20031220172052/http://paup.csit.fsu.edu/),
is invaluable in constructing these. Finally, write \"end;\" at the end
of the document. Save the document. Now, executing the document, the
same way you execute a data file, will start the commands running.

All the batch files on this website can be copied into a new PAUP
document and executed. Several of the batch files have places to insert
the name of a file or set, often appearing as [name]. Make sure that
the brackets are removed: for example, if the file is called tree1,
\"file=[name]\" should become \"file=tree1\". PAUP makes batch file
troubleshooting easy \-- if there is a problem with the file, PAUP will
open the file and put the cursor where the error occured. PAUP does not
see anything in brackets, so comments written there will not affect
analyses. The only two exceptions to this are: 1) if the first character
within the bracket is an exclamation point, PAUP will display the
contents in the display buffer; and 2) Trees can have a [&U] or [&R]
to tell PAUP if they are rooted or not.

[]{#ild}Testing for homogeneity, what to do if this fails

We use the partition homogeneity/incongruence-length difference test
implemented in PAUP to determine if different partitions of the data
(typically different loci) have significantly different signals. To do
this, you should first have a PAUP data file with a [character
partition](#charpartition) defined. Then, choose Partition Homogeneity
Test under the Analysis window to perform the test. Alternatively, you
could use the following batch file, which will start a log and perform
the test. You should replace [name] with the name of your partition
(remove the brackets).

> \#nexus\
> Begin PAUP;\
> log file=ildtest.log;\
> hompart partition=[name] nreps=100 / start=stepwise addseq=random
> nreps=10 savereps=no randomize=addseq rstatus=no hold=1 swap=tbr
> multrees=yes;\
> log stop;\
> end ;

This batch file sets the number of ILD replicates at 100, with the
number of random taxon addition replicates at 10 per ILD replicate. As
in everything else on this page, these are just starting suggestions,
not necessarily the best numbers for your data.

If the test is non-significant, there is no significant conflict between
the partitions. We then include all partitions in subsequent analyses.

If the test is significant, there may be conflict due to different
signals (horizontal transfer at one locus, for example) or due to a more
mundane problem in the data. To rule the latter out, first create all
possible charpartitions to compare pairs of loci. Then do the ILD test
again for each charparitition, excluding charsets not in the
charpartition being examined. For example, if we have three charsets
called 28s, coI, and ef1, a batch file would be:

> \#nexus\
> begin paup;\
> log file=pairwiseILD.log;\
> charpartition noef1 = 28s:28s, coI:coI;\
> charpartition no28s = coI:coI, ef1:ef1;\
> charpartition nocoI = 28s:28s, ef1:ef1;\
> exclude ef1 /only;\
> [!ILD between 28S and COI]\
> hompart partition=noef1 nreps=100;\
> exclude 28s /only;\
> [!ILD between COI and EF1]\
> hompart partition=no28s nreps=100;\
> exclude coI /only;\
> [!ILD between 28S and EF1]\
> hompart partition=nocoI nreps=100;\
> include all;\
> log stop;\
> end;

Note that the comments in brackets which start with an exclamation point
will be displayed when PAUP gets to that part of the batch file.

If comparisons between ef1 and either of the other two loci are
significant, but the comparison between 28sand coI are not significantly
different, you can begin to suspect that ef1 is different from the other
two genes. You could then go back to make sure that the sequences are
all labeled correctly as to taxon or that the sequences all come from
the same copy of the gene.

If a problem cannot be found in the sequences, the apparent conflict may
be real. Our lab\'s general approach is to analyze all the sequences
together, as we normally do, and do additional, separate analyses for
each locus. However, there is currently great debate within
phylogenetics on what to do in this case. Some people will analyze all
the data together regardless of apparent conflict, others will only
analyze together if there is no conflict.\

[]{#saturation}Saturation

Saturated data refers to data where the phylogenetic signal is
overwhelmed by multiple changes at each site. Note that data which is
saturated at deep levels may still be useful at more recent divergences.
Saturated data can be relatively harmless, merely making tree searches
take longer and potentially reducing the decisiveness of the dataset.
Some workers include all data, regardless of saturation, in the hopes
that there is still some phylogenetic signal amidst the noise. However,
saturated data can pose problems. The data could increase variance in
branch lengths, potentially raising the problem of long branch
attraction. This could also make clocklike data reject a clock. The
saturated data could require a more complex likelihood model than the
dataset would without the saturated data. Finally, the saturated data
could adversely affect the liklihood parameter estimates applied to all
sites, which may affect the tree chosen. For these reasons, we often
either exclude saturated data (i.e., excluding mitochondrial DNA third
codon positions) or convert it to less-saturated forms (i.e., converting
mitochondrial DNA sequence to amino acid sequence). Detecting saturated
data is a somewhat subjective process. One common approach is to graph
the transition/transversion ratio for pairs of sequences versus number
of transversions between those pairs of sequences. This will eventually
decline to an equilibrium value, governed by the base frequencies, at
saturation. Graphs for character sets which show a sharp decline and
long horizontal trend of data suggest saturation much more than graphs
which show a gradual decline in this ratio over all transversions.
Character sets with lower effective number of states (for example, third
codon positions of mitochondrial DNA in insects, which often shows
strong AT bias) are more likely to show saturation than character sets
with more balanced base composition.

[]{#exhaustive}Exhaustive and branch-and-bound searches (\<13 taxa)

These strategies are guaranteed to find the optimal tree but take too
long for more than a few taxa. As a result, we seldom have occassion to
use them.

[]{#heuristic}Heuristic searches

Heuristic searches are not guaranteed to find the globally-optimal tree,
but they can work for many more taxa than exhaustive or branch-and-bound
searches. A good starting [batch file](#batch) is:

> \#nexus\
> begin PAUP;\
> log file=hsearch1.log;\
> set autoclose=yes;\
> hsearch start=stepwise addseq=random nreps=**100** savereps=yes
> randomize=addseq rstatus=yes hold=1 swap=tbr multrees=yes;\
> savetrees file=hsearch1.all.tre;\
> filter best=yes permdel=yes;\
> savetrees file=hsearch1.best.tre;\
> log stop;\
> end;

This batch file will do a heuristic search, save all the trees found in
each random addition sequence replicate in the hsearch1.all.tre file,
then filter for the best trees overall and save them in an
hsearch1.best.tre file. PAUP will also output a tree-island profile. A
\"tree-island\" is a group of trees which cannot be reached through
branchswapping from a different group of trees. Ideally, all the trees
will be on one island which will have been hit every time, no matter
what starting tree or taxon addition order. If this is the case, the
treespace is fairly simple and we would probably move on to
[bootstrapping](#bootstrap). If the island(s) with the best trees was
hit nearly all the time, we\'d probably be satisfied, as well. However,
if the best trees are not recovered in many of the searches, we would
need to search further.

We use a variety of strategies for these next searches. The first
strategy is to simply increase the number of addition sequence
replicates. To do this, change the nreps=100 in the batch file above to
a higher number, perhaps nreps=1000. To avoid confusion, we may also
want to replace \"hsearch1\" with \"hsearch2\" wherever it appears
before running the search again.

Another strategy is to start from random trees, rather than random taxon
addition. On the above batch file, change the randomize=addseq to
randomize=trees, change the \"hsearch\#\" as above, and search again,
with nreps to taste.

These searches will run longer than the initial search. A way to speed
up the search while covering more of tree space is to increase the
number of addition replicates but reducing the thoroughness of the
search for each replicate. We might do this if we believe there may be
islands of good trees PAUP is missing. A sample batch file follows. You
may want to change the nreps and timelimit.

> \#nexus\
> begin paup;\
> log file=hsearch.tlimit.log;\
> set maxtrees=10000 increase=auto;\
> hsearch rstatus=no limitperrep=yes nreps=5000 randomize=trees
> timelimit=5 savereps=yes;\
> savetrees file=hsearch.tlimit.all.tre;\
> filter best=yes permdel=yes;\
> savetrees file=hsearch.tlimit.best.tre;\
> end;

The parsimony ratchet (below) can be useful in searching treespace
broadly, as well.


[]{#ratchet}Parsimony ratchet

The parsimony ratchet is a way to search treespace by reweights
characters for some iterations of a search. It is especially good for
searches with large numbers of taxa. It is described by Kevin Nixon
(Nixon, K. C. 1999. \"The Parsimony Ratchet, a new method for rapid
parsimony analysis.\" *Cladistics* 15: 407-414). Derek Sikes and Paul
Lewis have written a program, called
[PaupRat](https://web.archive.org/web/20031220172052/http://viceroy.eeb.uconn.edu/paupratweb/pauprat.htm),
which generates batch files to implement this search strategy in PAUP. A
few notes on the use of this: 1) It\'s better to do several searches of
a moderate number of nreps in each search (create a new ratchet file for
each search) than one search with many nreps; 2) It\'s useful to insert
the commands:

> stopcmd \"filter best=yes permdel=yes\";\
> stopcmd \"savetrees file=mydata.best.tre\";

into the setup.nex file before the stopcmd `quit`; to
automatically filter for the best trees.

[]{#bootstrap}Bootstrap searches

Bootstrap searches can take quite some time. One useful feature of PAUP
is the ability to search on multiple computers or on several different
occasions and combine the results. The key to this is saving the
bootstrap trees from each search and then loading them all together,
then computing the consensus tree using the tree weights. You must be
sure to keep the bootstrap search settings (except for the number of
bootstrap replicates) the same between searches for this to be valid.
This can be done through the menus (make sure to hit the \"save trees to
file\" checkbox in the bootstrap menu), or through the batch files
(below).

For the search itself he search can be stopped before completing all
the bootstrap replicates; if doing multiple searches, the treefile name
should be changed for each search. You may want to change the number of
bootstrap replicates (currently set at 500) and the number of random
taxon additions per bootstrap replicate (currently set at a low value of
10):

> \#nexus\
> begin paup;\
> set storetreewts=yes;\
> bootstrap nreps=**500** treefile=bootstrap1.tre replace=no/
> start=stepwise addseq=random nreps=**10** savereps=no randomize=addseq
> hold=1 swap=tbr multrees=yes;\
> end;

Load all the bootstrap trees, making sure to store tree weights (an
option which should have been made the default by the above batch file),
to load all blocks, and NOT to eliminate duplicate trees. Executing the
following batch file should set all these options and load the trees
saved as bootstrap1.tre in the active folder.

> \#nexus\
> begin paup;\
> gettrees allblocks=yes duptrees=keep storetreewts=yes mode=7
> file=bootstrap1.tre;\
> end;

Finally, compute a majority-rule consensus tree. The batch file for this
is:

> \#nexus\
> begin paup;\
> log file=bootstrapconsensus.log;\
