# Formulas to represent band edges in Kikuchi pattern

Herein, the origin is considered as the projection center, and the axes and the unit of the length is fixed so that the phosphor screen is included in the plane z=1. The following vectors are used in the following discussion:

- $`\mathbf{z}={(0,0,1)}^T`$
- a reciprocal lattice vector $`\mathbf{a}^*`$
- $`\mathbf{u}=\mathbf{a}^*/|\mathbf{a}^*|`$

The band center line in a Kikuchi pattern is the intersection of the phosphor screen $`\mathbf{z}\cdot{\mathbf{x}}=1`$ and the diffracting plane $`\mathbf{a}^*\cdot{\mathbf{x}}=0`$. The Bragg equation $`2d\sin{θ}=nλ`$ $`(n=1,2,...)`$ provides the relation between the Bragg angle θ and the d-spacing $`d=1/|\mathbf{a}^∗|`$.

For any point $`\mathbf{x}`$ of the Kossel Cone as in the following figure, $`(\mathbf{x}\cdot{\mathbf{u}})\mathbf{u}`$ is the projection of $`\mathbf{x}`$ in the direction of $`\mathbf{a}^*`$, and $`\mathbf{x}-(\mathbf{x}\cdot{\mathbf{u}})\mathbf{u}`$ is the projection in the plane orthogonal to $`\mathbf{a}^*`$. Hence, their lengths satisfy the below equality:

$`\dfrac{|(\mathbf{u}\cdot{\mathbf{x}})\mathbf{u}|^2}{|\mathbf{x}-(\mathbf{x}\cdot{\mathbf{u}})\mathbf{u}|^2}=\tan^2θ`$

![KosselCone2](https://github.com/arimuratak/ProjectEBSDConograph/assets/149344913/0c5e4dcf-5eaf-495b-b604-b21700cf5080)

Therefore, the Kossel cone is a conical surface defined by:

$`(\mathbf{u}\cdot{\mathbf{x}})^2=|\mathbf{x}|^2\sin^2θ`$

The band edges are the intersections of the Kossel cone and $`\mathbf{z}\cdot{\mathbf{x}}=1`$. For simplicity, we rotate the coordinate axes of the phosphor screen, and assume $`\mathbf{u}=(-\cos{σ},0,\sin{σ})`$. In this case, from $`\mathbf{u}\cdot{\mathbf{x}}=0`$, the equation of the band center line is provided by:

$`x=\tan{σ}`$

As for the equation of the band edges, we have:

$`(−x\cos{σ}+\sin{σ})^2=(x^2+y^2+1)\sin^2{θ}`$

As a result, the band edges are represented as the following hypebola:

$`(\cos^2{σ}-\sin^2{θ})(x-\dfrac{\cos{σ}\sin{σ}}{\cos^2{σ}-\sin^2{θ}})^2-(\sin^2{θ})y^2=\dfrac{\cos^2{θ}\sin^2{θ}}{\cos^2{σ}-\sin^2{θ}}`$

The above is also equivalent to:

$`(\dfrac{\cos{2σ}+\cos{2θ}}{\sin{2θ}})^2(x-\dfrac{\sin{2θ}}{\cos{2σ}+\cos{2θ}})^2-\dfrac{\cos{2σ}+\cos{2θ}}{2\cos^2{θ}}y^2=1`$

In the case of EBSD images, $`\cos{2σ}+\cos{2θ}>0`$ holds. The band width equals the distance between the two intersections $`σ_{begin}, σ_{end}`$ of the $`x`$-axis and the band edges: Furthermore, we have:

$`\tan{σ_{begin}}=\dfrac{\sin{2σ}-\sin{2θ}}{\cos{2σ}+\cos{2θ}}=
\dfrac{(e^{i(σ+θ)}+e^{-i(σ+θ)})(e^{i(σ-θ)}-e^{-i(σ-θ)})}{i(e^{i(σ+θ)}+e^{-i(σ+θ)})(e^{i(σ-θ)}+e^{-i(σ-θ)})}=\tan{(σ
-θ)}`$

$`\tan{σ_{end}}=\dfrac{\sin{2σ}+\sin{2θ}}{\cos{2σ}+\cos{2θ}}=
\dfrac{(e^{i(σ+θ)}-e^{-i(σ+θ)})(e^{i(σ-θ)}+e^{-i(σ-θ)})}{i(e^{i(σ+θ)}+e^{-i(σ+θ)})(e^{i(σ-θ)}+e^{-i(σ-θ)})}=\tan{(σ
+θ)}`$

As a result, $`σ_{begin}=σ−θ`$ , $`σ_{end}=σ+θ`$ are also obtained.
