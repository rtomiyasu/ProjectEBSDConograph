[to Japanese](https://github.com/rtomiyasu/ProjectEBSDConograph/blob/main/README-jp.md)
# Instructions for EBSD-CONOGRAPH (CUI version)

## Overview
This page explains the open source software [EBSD-CONOGRAPH Version 0.99](https://github.com/rtomiyasu/ProjectEBSDConograph/tree/main/EBSDConograph_0_9_99_win). The software executes ab-initio indexing of Kikuchi patterns (Figure 1) obtained by electron Backscatter Diffraction (EBSD). The indexing method is based on the CONOGRAPH method originally invented for [the powder auto-indexing software](https://github.com/rtomiyasu/ProjectPowderConograph).

![KikuchiPattern](https://github.com/rtomiyasu/ProjectEBSDConograph/assets/149344913/79144fc3-949f-4cda-84c1-c193fe564090)
```
Figure 1 : Kikuchi Pattern (Simulation data of steel) and a result of band extraction.
Yellow lines indicate band centers (or lines parallel to them).
The red segments, which are parts of the perpendicular lines through the pattern center,
indicate the band widths.
```
## Prior Information
Before starting the software, it is necessary to extract the following information from the Kikuchi pattern and input them from data.txt:
1. band-center positions (φ, σ),
1. band widths ($`σ_{begin}`$, $`σ_{end}`$).

EBSD indexing programs should also include codes for automatic PC calibration and band detection, but these functions have not been yet implemented into the distributed software. This also reflects the current situation that methods for each stage are still begin studied in order to improve their accuracy and reliability.

Methods and programs for ab-initio indexing can be developed independently from the above preprocessing stages. Even if the PC and band positions include errors to some degree, it is possible to estimate the unitcell parameters, although the following principle problems cannot be avoided:

- errors in the projection center, band center positions and band widths increase the errors and the uncertainties in the obtained unitcell parameters.
- there is an ambiguity in the ab-initio indexing (i.e., multiple distinct solutions can be as good as the correct one), which happens, in particular, in case of unclear band edges and unitcells with low symmetry [1].

## FAQ
- [How to use the EBSD CONOGRAPH program](#How_to_use_the_EBSD_CONOGRAPH_program)
- [On the contents of the input and output files](#On_the_contents_of_the_input_and_output_files)
  - [Input files](#Input_files)
  - [Output parameters](#Output_parameters)
  - [Figure of merit Mnew](#Figure_of_merit_Mnew)
- [When results are not satisfactory, which parameters in the input.txt should be modified?](#How_to_modify_parameters_in_input_txt)
- [How do I repot bugs?](#How_do_I_report_bugs)
- [How do I cite EBSD-CONOGRAPH?](#How_do_I_cite_EBSD_CONOGRAPH)

## How_to_use_the_EBSD_CONOGRAPH_program
1. The software EBSD-CONOGRAPH requires the following data.txt and input.txt as input files. (Examples can be found in the Sample folder.) 
    - input.txt: includes input parameters to adjust the search method and the output. ([Example](https://github.com/rtomiyasu/ProjectEBSDConograph/blob/main/EBSDConograph_0_9_99_win/sample/Fe(four_columns%2Cuse_only_band_centers)/input.txt))
    - data.txt: includes information about the band center positions and band widths.
        - Example 1: [data.txt](https://github.com/rtomiyasu/ProjectEBSDConograph/blob/main/EBSDConograph_0_9_99_win/sample/Fe(three_columns%2Cuse_band_widths)/data.txt) (3 columns: $`σ`$, $`σ_{begin}`$, $`σ_{end}`$. The σ is set to $`(σ_{begin} + σ_{end}) / 2`$) in this case,
        - Example 2: [data.txt](https://github.com/rtomiyasu/ProjectEBSDConograph/blob/main/EBSDConograph_0_9_99_win/sample/Fe3C(four_columns%2Cuse_band_width)/data.txt) (4 columns: $`φ`$, $`σ`$, $`σ_{begin}`$, $`σ_{end}`$).
        - The flag 0/1 in the first line of data.txt indicates either of the following options:
          - 1: estimation of the unitcell parameters from $`φ`$, $`σ`$, $`σ_{begin}`$, $`σ_{end}`$ <br>* The uncertainties of the unitcell parameters are increased in this case, because of errors in band widths (i.e., $`σ_{begin}`$, $`σ_{end}`$.)
          - 0: estimation of the length ratios a/c, b/c and the angles α, β, γ based on only φ and σ.
1. Copy one of the folders from the Sample folder. Modify the contents of data.txt (and input.txt in order to improve the indexing results).
1. Open a command prompt or terminal window in your operating system.<br>Change the current folder to the same folder that contains the modified input files.
1. Enter the path to the EBSDConograph.exe file on the Command Prompt window and execute EBSDConograph.

## On_the_contents_of_the_input_and_output_files
### Input_files
Figure 2 explains how $`φ`$, $`σ`$, $`σ_{begin}`$, $`σ_{end}`$ are provided by the band positions and widths in an EBSD image.
![KosselCone](https://github.com/rtomiyasu/ProjectEBSDConograph/assets/149344913/d944fc7c-c291-414b-830f-b5768005fba1)

- Figure 2:
  - (a) Band center lines can be regarded as the intersections of the phosphor screen and diffracting planes
through the projection center (PC).
The pattern center O is the foot of the perpendicular line from the PC to the phosphor screen.
  - (b) Band edges are the intersections of the phosphor screen and conical surfaces (the so called, Kossel cones).
They are parts of hyperbola lines ([formula](https://github.com/rtomiyasu/ProjectEBSDConograph/blob/main/html/FormulasForEBSDBandEdges_en.md)).
  - (c) The phosphor screen is parallel to the sheet.
The unit of the length is fixed so that the camera length (= distance between the PC and the screen) equals 1.
The φ is the angle between the x-axis and the perpendicular line from O to the band center line.
The $`σ`$, $`σ_{begin}`$, $`σ_{end}`$ are obtained from the distances of the center lines or the band edges from O.

The band widths are necessary in order to uniquely determine the ratios a/c, b/c and the angles α, β, γ, and obtain the unitcell scale (and hence, a, b, c). In fact, the Bragg angle θ is obtained from the band width by the calculation $`2θ = σ_{end} - σ_{begin}`$.

Futhermore, by the Bragg's law,

$`2d\sin{θ} = nλ`$, $`d`$ : spacing of the diffracting plane, $`n`$ : integer, $`λ`$ : wave length of the electron beam

Thus, if $`na^*`$($`n`$ : integer) are the reciprocal lattice vectors orthogonal to the diffracting plain,

$`\sin{θ} = nλ/2d = |na^*|λ/2`$, | | is the vector length.&emsp;&emsp;&emsp;(1)

As a result of Eq.(1), the bands correponding to the Miller indices $`n(hkℓ)`$ of $`na^∗`$ ($`n`$: integer), have the identical center lines but distinct band widths. Normally the band edges with $`n=±1`$ are the most clearly observed. However, visible bands are influenced by the magnitudes of structure factors and reflection rules due to the space-group symmetry (Nolze & Winkelmann, 2017), and furthermore, there seem to be exceptions (Fig. of Day (2008)).

### Output_parameters
The following are the output files of the software:
- Example 1: [out.txt](https://github.com/rtomiyasu/ProjectEBSDConograph/tree/main/EBSDConograph_0_9_99_win/sample/Fe(three_columns%2Cuse_band_widths)) (when band widths are also used),
- Example 2: [out.txt](https://github.com/rtomiyasu/ProjectEBSDConograph/blob/main/EBSDConograph_0_9_99_win/sample/Fe(four_columns%2Cuse_only_band_centers)/output/out.txt) (when only φ and σ are used).

In the output files, candidate unitcells are classified by their Bravais types, and sorted by the values of the figure of merit ($`M^{new}`$) defined in [1]. The following parameters are also output, after they are refined by using a non-linear squares method:

1. The basis $`a_1, a_2, a_3`$ of the crystal lattice that satisfies the following equation, and the estimated errors of their components:
   $` 
 \begin{pmatrix}
  a_1・a_1 & a_1・a_2 & a_1・a_3 \\
  a_2・a_1 & a_2・a_2 & a_2・a_3 \\
  a_3・a_1 & a_3・a_2 & a_3・a_3 
 \end{pmatrix} = \begin{pmatrix}
  a^2 & ab\cos{γ} & ac\cos{β} \\
  ab\cos{γ} & b^2 & bc\cos{α} \\
  ac\cos{β} & bc\cos{α} & c^2 
 \end{pmatrix}, a,b,c,α,β,γ`$ : unitcell parameters.

2. Euler angles $`θ_1`$, $`θ_2`$, $`θ_3`$ that represents the direction of the lattice (more precisely, the following orthogonal matrix $`G`$
), and their estimated errors:<br>
   $`G:=L^{-1}A=
   \begin{pmatrix}\cos{θ_1} & \sin{θ_1} & 0 \\-\sin{θ_1} & \cos{θ_1} & 0 \\0 & 0 & 1\end{pmatrix}
   \begin{pmatrix}1 & 0 & 0 \\0 & \cos{θ_2} & \sin{θ_2} \\0 & -\sin{θ_2} & \cos{θ_2} \end{pmatrix}
   \begin{pmatrix}\cos{θ_3} & \sin{θ_3} & 0 \\-\sin{θ_3} & \cos{θ_3} & 0 \\0 & 0 & 1\end{pmatrix}`$,

   where $`A`$ is the 3×3 matrix with the row vectors $`a_1`$, $`a_2`$, $`a_3`$, and $`L`$ is the is the lower triangle matrix with $`LL^T = AA^T`$.

3. Projection-center shift $`Δx`$, $`Δy`$, $`Δz`$, and their estimated errors. If the camera length and the pattern-center coordinate used to make the input file data.txt are denoted by Lobs
 and ($`X^{old}`$, $`Y^{old}`$), the new camera length $`L^{new}`$ and the pattern-center $`(X^{new}, Y^{new})`$ coordinate are given by:

   $`L^{new}=(1-Δz)L^{old}`$
   
   $`(X^{new}, Y^{new}) = (X^{old}, Y^{old}) + (L^{old}Δx, L^{old}Δy)`$

The following should be noted, with regard to the above parameters:
- in some cases, the parameters become worse by the refinement process, as a result of changes in the Miller indices assigned to the bands.
- the above estimated errors are propagation errors, when it is assumed in the non-linear least squares method that the input angles $`φ`$, $`σ`$, $`σ_{begin}`$, $`σ_{end}`$ have errors within 1 degree.

### Figure_of_merit_Mnew
Since $`M^{new}`$ is a generalization of the de Wolff M used in powder indexing, $`M`$ and $`M^{new}`$ have similar properties, except that M prefers higher-symmetric cells ([c.f.](https://github.com/rtomiyasu/ProjectEBSDConograph/blob/main/figures/table5_2_en.png)). Namely,
- If $`M^{new} > 10`$ for some unitcells, there are much possibility that the correct solution has been obtained.
- The correct solution should be one of the unitcells that gained approximately the largest $`M^{new}`$ value.
- Similar unitcell parameters have similar $`M^{new}`$ values, even if they belong to distinct Bravais types.

Therefore, user should check both the value of the figure of merit and the Bravais types, when selecting the correct solution from the output list.

## How_to_modify_parameters_in_input_txt
The following two options can be chosen as a searching method:
- Quick search (often enough for unitcells with high symmetry),
- Exhaustive search (for all cases).

The searched region can be also expanded by increasing the following parameters (the upper ones are more influential. Also see explanations in input.txt), although comparatively large numbers are already used in the attached input.txt:

- Upper bound on errors in phi, sigma, sigma_begin, sigma_end: 1 degree
- Max |h|,|k|,|l| used for indexing : 7
- Tolerance level for errors in the unitcell scales : 3
- Resolution for Bravais-type determination : 0.02

## How_do_I_report_bugs
You should send us a bug report with all the input and output files attached (including LOG_CONOGRAPH.txt) to the following e-mail address:

- tomiyasu.ryoko.446 (at) m.kyushu-u.ac.jp

## How_do_I_cite_EBSD_CONOGRAPH
If you use the program for your research, we strongly encourage you to include a citation of the following article in the bibliography.

- R. Oishi-Tomiyasu, T. Tanaka, J. Nakagawa, “Distribution rules of systematic absence and generalized deWolff figure of merit applied to EBSD ab-initio indexing”, [arxiv](https://arxiv.org/abs/2003.13403).
