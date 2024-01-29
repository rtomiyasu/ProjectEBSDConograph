[to_English](https://github.com/rtomiyasu/ProjectEBSDConograph/blob/main/html/FormulasForEBSDBandEdges_en.md)

# 菊池パターンのバンドエッジを表す式

原点をプロジェクションセンターとし、蛍光板がz=1の平面内にくるように、長さの単位を設定する。 以下のベクトルを用いて式を導出する:

- $`\mathbf{z}={(0,0,1)}^T`$
- 逆格子ベクトル $`\mathbf{a}^*`$
- $`\mathbf{u}=\mathbf{a}^*/|\mathbf{a}^*|`$

バンドセンターは、 回折面に対応する平面 $`\mathbf{a}^*\cdot{\mathbf{x}}=0`$ と 蛍光版 $`\mathbf{z}\cdot{\mathbf{x}}=1`$ の交わりになる。ブラッグの式 $`2d\sin{θ}=nλ`$ $`(n=1,2,...)`$ により、ブラッグ各θは、格子面間隔 $`d=1/|\mathbf{a}^∗|`$ と関係づけられる。

以下の図のKossel coneに含まれる座標 $`\mathbf{x}`$ は、 その $`\mathbf{a}^*`$ 方向への射影 $`(\mathbf{x}\cdot{\mathbf{u}})\mathbf{u}`$ とその鉛直方向への射影 $`\mathbf{x}-(\mathbf{x}\cdot{\mathbf{u}})\mathbf{u}`$ を考えたとき、下記の等式が成立する。

$`\dfrac{|(\mathbf{u}\cdot{\mathbf{x}})\mathbf{u}|^2}{|\mathbf{x}-(\mathbf{x}\cdot{\mathbf{u}})\mathbf{u}|^2}=\tan^2θ`$

![KosselCone2](https://github.com/rtomiyasu/ProjectEBSDConograph/assets/149344913/bca08665-6dce-42b2-ac37-5acd52d4fc7f)

したがって、Kossel coneは以下で定義される円錐面になる:

$`(\mathbf{u}\cdot{\mathbf{x}})^2=|\mathbf{x}|^2\sin^2θ`$

バンドエッジはこのKossel coneと平面 $`\mathbf{z}\cdot{\mathbf{x}}=1`$ の交わりになる。 簡単のため、蛍光板の座標軸を回転して、$`\mathbf{u}=(-\cos{σ},0,\sin{σ})`$ としておく。 $`\mathbf{u}\cdot{\mathbf{x}}=0`$ より、バンドセンターの方程式は以下に等しい。

$`x=\tan{σ}`$

バンドエッジの方程式は、

$`(−x\cos{σ}+\sin{σ})^2=(x^2+y^2+1)\sin^2{θ}`$

すなわち、以下の双曲線になる。

$`(\cos^2{σ}-\sin^2{θ})(x-\dfrac{\cos{σ}\sin{σ}}{\cos^2{σ}-\sin^2{θ}})^2-(\sin^2{θ})y^2=\dfrac{\cos^2{θ}\sin^2{θ}}{\cos^2{σ}-\sin^2{θ}}`$

さらに整理すると、

$`(\dfrac{\cos{2σ}+\cos{2θ}}{\sin{2θ}})^2(x-\dfrac{\sin{2θ}}{\cos{2σ}+\cos{2θ}})^2-\dfrac{\cos{2σ}+\cos{2θ}}{2\cos^2{θ}}y^2=1`$

今、$`\cos{2σ}+\cos{2θ}>0`$ は仮定してよい。バンドエッジと$`x`$軸の、2つの交点の距離 $`σ_{begin}, σ_{end}`$ がバンド幅に該当する。 さらに以下が成立する。 

$`\tan{σ_{begin}}=\dfrac{\sin{2σ}-\sin{2θ}}{\cos{2σ}+\cos{2θ}}=
\dfrac{(e^{i(σ+θ)}+e^{-i(σ+θ)})(e^{i(σ-θ)}-e^{-i(σ-θ)})}{i(e^{i(σ+θ)}+e^{-i(σ+θ)})(e^{i(σ-θ)}+e^{-i(σ-θ)})}=\tan{(σ
-θ)}`$

$`\tan{σ_{end}}=\dfrac{\sin{2σ}+\sin{2θ}}{\cos{2σ}+\cos{2θ}}=
\dfrac{(e^{i(σ+θ)}-e^{-i(σ+θ)})(e^{i(σ-θ)}+e^{-i(σ-θ)})}{i(e^{i(σ+θ)}+e^{-i(σ+θ)})(e^{i(σ-θ)}+e^{-i(σ-θ)})}=\tan{(σ
+θ)}`$

したがって、$`σ_{begin}=σ−θ`$ , $`σ_{end}=σ+θ`$ が得られる。
