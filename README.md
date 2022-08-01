# lstmContacts
The **lstmContacts** software characterizes and draws atigen-antibody interaction profiles based on a library of known contacts. The software has two main modules: the contact similarity search module (R module) and the prediction module, based on long-short term memory (LSTM) recurrent neural networks (Python module). The R module takes the input antigen-antibody complex and generates an affinity contact library that can be used by the LSTM module for the refinement/prediction of an interaction time series.

## 1.1. Requirements and installation

...

## 1.2. Internal representation of an antigen-antibody complex

The input antigen-antibody **complex** (AAC) is represented by a list x = [x1, x2, …, xn] of n **contacts**. Each j-th contact is a vector xj = [x0j, x1j, …, xmj] of mj elements, where x0j is the antibody residue interacting with the x1j, …, xmj residues on the surface of the receptor binding domain (RBD) of the Spike protein variant. Every contact corresponds univocally to a vector ai = affinity(xj) of 101 affinity score values, ranging from 0 to 1, and corresponding to the 101 nanoseconds of the molecular dynamics simulation stored in the internal library (object `contact.data`). In the internal library, antibodies are reported with the corresponding protein data bank ([**PDB**](https://www.rcsb.org/)) 3D structure ID: 7kmg (Bamlanivimab, Ly-Cov555), 7c01 (Etesevimab, Ly-Cov016), 7cm4 (Regdanvimab, CTP-59), 7l7d (Tixagevimab, AZD8895), 7l7e (Cilgavimab, AZD1061), 7r6w (Sotrovimab), 6zcz (EY6A). The variants of the library for which a molecular dynamic simulation is available include: *wt*, *alpha*, *beta*, *delta*, *omicron*. The input AAC wil be searched, both exactly and by similarity, against these data.

Computationally, the AAC must be specified as follows (R code):

```r
x <- list(x1 = c("h.R50", "V483", "E484"),
          x2 = c("h.L55", "L452", "T470", "F490"),
          x3 = c("h.Y101", "E484", "F490"),
          x4 = c("h.R104", "Q493", "S494"),
          x5 = c("l.Y32", "F486", "Y489"),
          x6 = c("l.Y92", "F486", "Y489"),
          x7 = c("l.R96", "V483", "E484"))
```

In the example above, `x` is the AAC and each vector in `x` is a contact. The first element of each contact is always the amino acid residue of the antibody, defined by the FAB chain ("h" for "heavy" and "l" for "light"), followed by a dot, the single letter code of the residue, and its position in the polypeptide chain. The other elements of a contact are the antigen residues interacting with the antibody one.

The `preprocess(vj)` function of the R module assigns the best possible ai vector to the j-th contact. If the dynamics for vj are in our contact library, the assignment is referred to as an *exact match*. Otherwise, the [x0j, x1j, …, xmj] residues of the vj contact are searched by similarity. Each residue xkj (with k = 0, 1, …, m) is firstly searched in a nearby position of the polypeptide chain. If this search fails, the algorithm seeks for a residue with similar chemical properties (referred to as "group"), according to the table below (`contact.groups` object). Since the dynamics are run in triplicates, the internal library has three affinity vector replicates per contact: ai1, ai2, and ai3, where ai is computed as the element-wise mean of the three vectors.

```
> contact.groups
               aa tlc code     group                         description
1         Alanine Ala    A aliphatic                    Nonpolar neutral
2      Isoleucine Ile    I aliphatic                    Nonpolar neutral
3         Leucine Leu    L aliphatic                    Nonpolar neutral
4          Valine Val    V aliphatic                    Nonpolar neutral
5      Methionine Met    M aliphatic Nonpolar neutral, sulfur-containing
6   Phenilalanine Phe    F  aromatic             Nonpolar neutral, large
7     Thryptophan Trp    W  aromatic             Nonpolar neutral, large
8        Tyrosine Tyr    Y  aromatic            Polar hydrophobic, large
9      Asparagine Asn    N     polar                     Polar uncharged
10      Glutamine Glu    Q     polar                     Polar uncharged
11         Serine Ser    S     polar                     Polar uncharged
12      Threonine Thr    T     polar                     Polar uncharged
13       Arginine Arg    R     basic                         Charged (+)
14      Histidine His    H     basic                         Charged (+)
15         Lysine Lys    K     basic                         Charged (+)
16      Aspartate Asp    D    acidic                         Charged (-)
17      Glutamate Glu    E    acidic                         Charged (-)
18        Glycine Gly    G      tiny          Hydrogen side chain (tiny)
19        Proline Pro    P    cyclic             Cyclic nonpolar neutral
20       Cysteine Cys    C   sbridge    Charged (-), SH-containing polar
21 Selenocysteine Sec    U  selenium   Charged (-), SeH-containing polar
22    Pyrrolysine Pyr    O     amber   Charged (+), Bacteria and Archaea
```
For the antibody residue x0j of vj, the similarity search is further constrained within the original FAB chain (either *heavy* or *light*). A warning level based on powers of 2 is used to determine the quality of the search results: 2^3 = 8, exact match failed; 2^2 = 4, residue match (with the same FAB chain as the input) failed; 2^1 = 2, group search failed; 2^0 = 1, residue match (different FAB chain from the input) failed. A warning level below 8 means that an exact match was found, whereas a warning level from 8 to 14 indicates a non-exact match. If the warning level reaches 15, the search failed at every level and the algorithm cannot go further. This warning system enables fast user monitoring of the prediction quality. By default, the R module shows the search status, commenting on results quality in human language.

Residue search is further refined by two criteria. Firstly, each contact is classified based on the available interaction class. The class of an interaction specifies which part of the surface is involved in the AAC formation. By default, the object `contact.class` specifies the four available classes (I-to-IV) for the Spike RBD-FAB interaction:

```
> contact.class
$I
class1 class1 class1 class1 class1 class1 class1 class1 class1 class1 class1 
"R403" "D405" "E406" "R408" "Q409" "T415" "G416" "K417" "D420" "Y421" "Y453" 
class1 class1 class1 class1 class1 class1 class1 class1 class1 class1 class1 
"L455" "F456" "R457" "K458" "N460" "Y473" "A475" "G476" "S477" "F486" "N487" 
class1 class1 class1 class1 class1 class1 class1 class1 class1 class1 
"Y489" "Q493" "S494" "Y495" "G496" "Q498" "T500" "N501" "G502" "Y505" 

$II
class2 class2 class2 class2 class2 class2 class2 class2 class2 class2 class2 
"R403" "D405" "K417" "G446" "Y449" "N450" "L452" "L455" "F456" "T470" "S477" 
class2 class2 class2 class2 class2 class2 class2 class2 class2 class2 class2 
"N481" "V483" "E484" "G485" "F486" "N487" "Y489" "F490" "L492" "Q493" "S494" 
class2 class2 
"Q498" "Y505" 

$III
class3 class3 class3 class3 class3 class3 class3 class3 class3 class3 class3 
"N334" "P337" "E340" "T345" "R346" "R357" "N440" "L441" "K444" "V445" "G446" 
class3 class3 class3 class3 class3 class3 class3 class3 class3 class3 class3 
"G447" "N448" "Y449" "N450" "L452" "R466" "T470" "E471" "E484" "F490" "Q493" 
class3 class3 class3 
"Q498" "P499" "Y505" 

$IV
class4 class4 class4 class4 class4 class4 class4 class4 class4 class4 class4 
"Y369" "N370" "S371" "A372" "S373" "F374" "S375" "T376" "F377" "K378" "C379" 
class4 class4 class4 class4 class4 class4 class4 class4 class4 class4 class4 
"Y380" "G381" "V382" "S383" "P384" "T385" "K386" "L390" "D405" "R408" "Q409" 
class4 class4 class4 class4 class4 class4 class4 class4 class4 class4 
"P412" "Q414" "D427" "D428" "F429" "T430" "N437" "V503" "G504" "Y508" 
```

Each residue is assigned (also by similarity) to one or more classes. Then the odds ratio p = odds(i-th class)/odds(not i-th class) is calculated. The contact is assigned to the class with the highest p, referred to as the "majority class". When all the contacts of a complex belong to the same class, we have a confirmation to the goodness of the search results. 
Secondly, to find the antigen variants that best fit the input complex, raw results are further refined by checking for possible mutation hotspots. The default object `contact.mutations` contains mutational information for the Spike variants alpha, beta, delta, and omicron, with respect to the wt (wild-type) protein:

```
> contact.mutations
   variant   wt mutant  wt.group mutant.group
1    alpha N501   Y501     polar     aromatic
2     beta K417   N417     basic        polar
3     beta E484   K484    acidic        basic
4     beta N501   Y501     polar     aromatic
5    delta L452   R452 aliphatic        basic
6    delta T478   K478     polar        basic
7  omicron G339   D339      tiny       acidic
8  omicron S371   L371     polar    aliphatic
9  omicron S373   P373     polar       cyclic
10 omicron S375   F375     polar     aromatic
11 omicron K417   N417     basic        polar
12 omicron N440   K440     polar        basic
13 omicron G446   S446      tiny        polar
14 omicron S477   N477     polar        polar
15 omicron T478   K478     polar        basic
16 omicron E484   A484    acidic    aliphatic
17 omicron Q493   R493     polar        basic
18 omicron G496   S496      tiny        polar
19 omicron Q498   R498     polar        basic
20 omicron N501   Y501     polar     aromatic
21 omicron Y505   H505  aromatic        basic
```

## 1.3. Contact search

First, we need to open an R console and load the needed functions (lstmContacts requires only the basic R environment, version >= 4.0):

```r
source("~/lstmContacts/contacts.R")
```

The lstmContacts software allows to search for a single contact, through the `preprocess(x)` function:

```r
# Define a contact
x <- c("h.R105", "Y449", "L455", "L492", "Q493", "S494")

# Search
R <- preprocess(x)
```

The `preprocess` function prints out a summary composed by four messages:

```
### p
1.09261696547824
4.36452855190691
0.870856650882479
0.007049945387747
### Majority class: 2
###  Warning level: 0/15 [ok]
###  Contact match: exact
```

- the value of the interaction class odds-ratio (p);
- the selected majority class based on p;
- the warning level (0-7 [ok], 8-14 [suboptimal], 15 [unreliable]);
- the type of search result (*exact*: perfect contact match found in the library, *similarity*: output contact defined by similarity).

The search results are stored as a list:

```
> R
$input
[1] "h.R105" "Y449"   "L455"   "L492"   "Q493"   "S494"  

$antibody
[1] "7cm4"

$ab.residues
    7cm4 
"h.R105" 

$exact.match
  7cm4   7cm4   7cm4   7cm4   7cm4 
"Y449" "L455" "L492" "Q493" "S494" 

$similarity.match
[1] "None"

$class.residues
[1] "Y449" "L455" "L492" "Q493" "S494"

$class
[1] "class2"

$warning.level
[1] 0

$contact.level
[1] 3

$variant
[1] "omicron"

$hotspot
[1] "None"
```

The search output list includes: a copy of the input contact (`input`), the candidate antibody (`antibody`), which of the residues of the chosen antibodies best match the input (`ab.residues`), the results of the exact match for the input antigen and interacting antibody (`exact.match`; "None" if no exact match is found), the results of the similarity match (`similarity.match`; disabled if an exact match was found), the residues of the internal library matching those of the antigen variant (`class.residues`), the majority class (`class`), the warning level of the search result (`warning.level`; from 0 to 15), quality of the search (`contact.level`; 3: exact, 2: good similarity search, 1: suboptimal similarity search), the proposed variant (`variant`), possible mutant residues (`hotspot`).

Let us make an example with a mock contact, to see how the similarity search output looks like:

```r
# Define a contact
x <- c("h.D105", "Q449", "W455", "M492", "L493", "K494")

# Search
R <- preprocess(x)
```

As we see, the summary raises some warning messages:

```
### p
0.628126310459145
0.925496248203927
2.423437389895
0.628126310459145
### Majority class: 3
###  Warning level: 11/15 [suboptimal]
###  Contact match: similarity
Messaggi di avvertimento:
1: In ab.search(x[1], ab.library) : Partial antibody match.
2: In res.search(res, x.ag) : Modeling possible suboptimal contacts.
```

If we inspect the results, ...

```
> R
$input
[1] "h.D105" "Q449"   "W455"   "M492"   "L493"   "K494"  

$antibody
[1] "7cm4" "7l7e" "6zcz" "7r6w"

$ab.residues
    7cm4     7cm4     7cm4     7cm4     7l7e     7l7e     6zcz     6zcz 
 "h.D54"  "h.D56"  "h.D57"  "h.D51"  "h.D56" "h.D107"  "h.D33"  "h.D99" 
    7r6w 
"h.E108" 

$similarity.match
  7cm4   7l7e   6zcz   7r6w 
"Q493" "K444" "K386" "K356" 

$class.residues
[1] "K444" "Q493"

$class
[1] "class3"

$warning.level
[1] 11

$contact.level
[1] 2

$variant
[1] "omicron" "wt"     

$hotspot
[1] "W455" "M492" "L493"

$exact.match
[1] "None"
```

1.4. AAC search

The user can search for an entire complex (AAC), specifying it as reported in section 1.2. Let us see a quick example:

