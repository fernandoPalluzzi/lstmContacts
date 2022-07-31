# lstmContacts
The **lstmContacts** software characterizes and draws atigen-antibody interaction profiles based on a library of known contacts.

## 1.1. Installation

...

## 1.2. Internal representation of an antigen-antibody complex

The **lstmContacts** software has two main modules: the contact similarity search module (R module) and the prediction module, based on long-short term memory (LSTM) recurrent neural networks (Python module). The R module takes the input antigen-antibody complex and generates an affinity contact library that can be used by the LSTM module for the prediction.

The input **complex** is represented by a list C = [v1, v2, …, vn] of n **contacts**. Each j-th contact is a vector vj = [x0j, x1j, …, xmj] of mj elements, where x0j is the antibody residue interacting with the x1j, …, xmj residues on the surface of the receptor binding domain (RBD) of the Spike protein variant. Every contact corresponds univocally to a vector ai = affinity(vj) of 101 affinity score values, ranging from 0 to 1, and corresponding to the 101 nanoseconds of the molecular dynamics simulation. Since the dynamics are run in triplicates, we have three affinity vector replicates per contact: ai1, ai2, and ai3, where ai is computed as the element-wise mean of the three vectors.
The preprocess(vj) function of the R module assigns the best possible ai vector to the vj contact. If the dynamics for vj are in our contact library, the assignment is referred to as an exact match. Otherwise, the [x0j, x1j, …, xmj] residues of the vj contact are searched by similarity. Each residue xkj (with k = 0, 1, …, m) is firstly searched in a nearby position of the polypeptide chain. If this search fails, the algorithm seeks for a residue with similar chemical properties (referred to as “group”), according to the table below (`contact.groups` object).

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
For the antibody residue x0j of vj, the similarity search is further constrained within the original FAB chain (either *heavy* or *light*). A warning level based on powers of 2 is used to determine the search results quality: 2^3 = 8, exact match failed; 2^2 = 4, residue match (same FAB chain as the input) failed; 2^1 = 2, group search failed; 2^0 = 1, residue match (FAB chain mismatch) failed. A warning level below 8 means that an exact match was found, whereas a warning level > 8 and < 15 indicates a non-exact match. If the warning level reaches 15, the search failed at every level and the algorithm cannot go further. This warning system enables fast user monitoring of the prediction quality. By default, the R module shows the search status, commenting on results quality in human language.

## 1.3. Contact search and complex definition

...
