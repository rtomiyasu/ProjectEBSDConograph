[to English](https://github.com/rtomiyasu/ProjectEBSDConograph/blob/main/README.md)

# EBSD-CONOGRAPH (CUI版プログラム)の操作説明書

## 概要
以下では，オープンソースの[EBSD-CONOGRAPH Version 0.99](https://github.com/rtomiyasu/ProjectEBSDConograph/tree/main/EBSDConograph_0_9_99_win) について説明します．
このプログラムは， もともと粉末回折のために開発された[Conograph](https://github.com/rtomiyasu/ProjectPowderConograph)の方法に基づき， 電子線後方回折(EBSD)の菊池パターン (Figure 1)のバンド情報から，格子定数の推定とバンドの指数付けを同時に行うab-initio indexingを実施します．

![KikuchiPattern](https://github.com/rtomiyasu/ProjectEBSDConograph/assets/149344913/79144fc3-949f-4cda-84c1-c193fe564090)
```
Figure 1 : 菊池パターン(鋼のシミュレーションデータ)とバンド抽出結果．黄色線はバンドの中心線（またはその平行線）を示す．
赤線はパターンセンターを通るバンドの垂線の一部で，バンド幅を示す．
```
## NEWS
### 2024/1/31
- 配布ページを OSDN から GitHub に移動しました。

## 事前情報
ソフトウェア実行前に，菊池パターンから抽出する必要がある情報は:
1. バンドの中心線の座標 (φ, σ),
1. バンド幅 $`(σ_{begin}`$, $`σ_{end})`$.

EBSDソフトウェアには，PC座標の補正とバンド検出を自動で行うプログラムが必要ですが， 配布中のソフトウェアにはまだ実装されていません． これは，それぞれの手法の精度および信頼性の向上が研究途上であることも反映しています．

指数付けの手法およびプログラムは，上記の事前処理とは独立に実装できます． PCやバンド座標がある程度の誤差を含む場合も，以下の原理的な問題を無視すれば， 正解に近い格子定数を得ることは可能です:

- PCやバンドの中心線，バンド幅の誤差が大きいほど，得られたユニットセルの精度も悪くなる.
- ab-initio indexingでは，同程度に良い複数の異なる解が生じることがあり， 特に，バンドエッジが明瞭でなく格子の対称性が低いときに起こる [1]．

## FAQ
- [EBSD-CONOGRAPHの使い方](#EBSD-CONOGRAPHの使い方)
- [入出力ファイルの内容について](#入出力ファイルの内容について)
  - [入力ファイル](#入力ファイル)
  - [出力されるパラメータ](#出力されるパラメータ)
  - [Figure of merit Mnew](#Figure_of_merit_Mnew)
- [うまくいかないときに変更するinput.txtのパラメータ](#うまくいかないときに変更するinput_txtのパラメータ)
- [バグ報告のお願い](#バグ報告のお願い)
- [論文の引用に関して](#論文の引用に関して)

## EBSD-CONOGRAPHの使い方
1. 同プログラムを実行するには，以下のdata.txt, input.txtを入力ファイルとして準備する必要があります． （付属の Sample フォルダに入力例があります．）
    - input.txt: 探索ファイルや出力を調整する入力パラータを含む ([例](https://github.com/rtomiyasu/ProjectEBSDConograph/blob/main/EBSDConograph_0_9_99_win/sample/Fe(four_columns%2Cuse_only_band_centers)/input.txt))。
    - data.txt: バンドの中心線とバンド幅の情報を含む．
        - Example 1: [data.txt](https://github.com/rtomiyasu/ProjectEBSDConograph/blob/main/EBSDConograph_0_9_99_win/sample/Fe(three_columns%2Cuse_band_widths)/data.txt) (3列データ: $`φ`$, $`σ_{begin}`$, $`σ_{end}`$。 このときσは、$`(σ_{begin} + σ_{end}) / 2`$とする),
        - Example 2: [data.txt](https://github.com/rtomiyasu/ProjectEBSDConograph/blob/main/EBSDConograph_0_9_99_win/sample/Fe3C(four_columns%2Cuse_band_width)/data.txt) (4列データ: $`φ`$, $`σ`$, $`σ_{begin}`$, $`σ_{end}`$).
        - data.txtの1行目では，以下のいずれを実行するかのフラグ 0/1を指定してください(3,4列データどちらでも同じ)．
          - 1: $`φ`$, $`σ`$, $`σ_{begin}`$, $`σ_{end}`$からの，格子定数の推定 * この場合，一般に，バンド幅($`σ_{begin}`$, $`σ_{end}`$)の精度が悪いため，得られる格子定数の精度も悪くなります。
          - 0: $`φ`$, $`σ`$からの，比 $`a/c`$, $`b/c`$ と角度 $`α`$, $`β`$, $`γ`$ の推定．
1. Sampleフォルダの中の一つのフォルダをコピーし， コピー先のフォルダ内のdata.txtを(解析結果をもっとよくしたい場合，input.txtも)適宜修正してください。
1. コマンドプロンプト，またはお使いのOSのターミナルウィンドウを立ち上げ， 上で変更した入力ファイルと同じフォルダにカレントフォルダを移動してください
1. コマンドプロンプトより，EBSDConograph.exeのパスを入力して実行します．

## 入出力ファイルの内容について
### 入力ファイル
Figure 2 にexplains how $`φ`$, $`σ`$, $`σ_{begin}`$, $`σ_{end}`$ がEBSD画像のバンド座標, バンド幅からどのように与えられるかを示します．
![KosselCone](https://github.com/rtomiyasu/ProjectEBSDConograph/assets/149344913/d944fc7c-c291-414b-830f-b5768005fba1)

- Figure 2:
  - (a) バンドの中心線は， 蛍光板とprojection center(PC)を通る回折面の交わり。パターンセンター O はPCから蛍光板に下した垂線の足の座標に等しい。
  - (b) バンドエッジは，蛍光板と円錐面 (Kossel cone)の交わり，よって双曲線になる ([式](https://github.com/rtomiyasu/ProjectEBSDConograph/blob/main/html/FormulasForEBSDBandEdges_jp.md)). 
  - (c) 蛍光板は紙面と並行とし，長さの単位はカメラ長(= PCと蛍光板の距離)が1となるよう設定する． このとき，φは，X軸とPC2からバンドの中心線に下した垂線がなす角度に等しい． また，Oとバンドの中心線, バンドエッジとの距離から $`σ`$, $`σ_{begin}`$, $`σ_{end}`$が得られる．

バンド幅情報は，長さの比a/c, b/cと角度α, β, γ を一意に決めるため，さらに格子定数のスケール (よってa, b, c)を得るために必要です．実際，図2(b)にみるように，ブラッグ角 θ は， $`2θ = σ_{end} - σ_{begin}`$ より，バンド幅から得ることができます．さらにブラッグの法則より，

$`2d\sin{θ} = nλ`$, $`d`$ : 回折面の格子面間隔, $`n`$ : 整数, $`λ`$ : 電子線ビームの波長

このとき，各回折面に垂直な逆格子ベクトルを $`na^*`$($`n`$ : 整数)とすると、（||はベクトルの長さ）

$`\sin{θ} = nλ/2d = |na^*|λ/2`$&emsp;&emsp;&emsp;(1)

式(1)より、$`na^∗`$ ($`n`$: 整数)のミラー指数 $`n(hkℓ)`$ に バンドの中心線は同一ですが， バンド幅は異なります。もっともよく見えるのは $`n=±1`$ のバンドであることが多いですが， 構造因子の大きさや空間群の消滅則の影響があり(Nolze & Winkelmann, 2017)，さらに例外もあるようです(Fig.19 of Day (2008))。

### 出力されるパラメータ
以下はソフトウェアの出力ファイルの例です:
- Example 1: [out.txt](https://github.com/rtomiyasu/ProjectEBSDConograph/blob/main/EBSDConograph_0_9_99_win/sample/Fe(three_columns%2Cuse_band_widths)/output/out.txt) (バンド幅を用いた場合)
- Example 2: [out.txt](https://github.com/rtomiyasu/ProjectEBSDConograph/blob/main/EBSDConograph_0_9_99_win/sample/Fe(four_columns%2Cuse_only_band_centers)/output/out.txt) (φ, σのみを用いた場合)

出力ファイルout.txtの各格子定数は，ブラべー格子で分類され， [1]で定義されたfigure of merit ($`M^{new}`$)でソートされています。 また以下のパラメータも非線形最小二乗法による精密化の後、出力されます。

1. 以下の関係式を満たす実格子基底 $`a_1, a_2, a_3`$ とその各成分の推定誤差::
   $` 
 \begin{pmatrix}
  a_1・a_1 & a_1・a_2 & a_1・a_3 \\
  a_2・a_1 & a_2・a_2 & a_2・a_3 \\
  a_3・a_1 & a_3・a_2 & a_3・a_3 
 \end{pmatrix} = \begin{pmatrix}
  a^2 & ab\cos{γ} & ac\cos{β} \\
  ab\cos{γ} & b^2 & bc\cos{α} \\
  ac\cos{β} & bc\cos{α} & c^2 
 \end{pmatrix}, a,b,c,α,β,γ`$ : 格子定数
2. 実格子の向き(具体的には以下の直交行列 $`G`$)を表すオイラー角 $`θ_1`$, $`θ_2`$, $`θ_3`$ とその推定誤差:<br>
   $`G:=L^{-1}A=
   \begin{pmatrix}\cos{θ_1} & \sin{θ_1} & 0 \\-\sin{θ_1} & \cos{θ_1} & 0 \\0 & 0 & 1\end{pmatrix}
   \begin{pmatrix}1 & 0 & 0 \\0 & \cos{θ_2} & \sin{θ_2} \\0 & -\sin{θ_2} & \cos{θ_2} \end{pmatrix}
   \begin{pmatrix}\cos{θ_3} & \sin{θ_3} & 0 \\-\sin{θ_3} & \cos{θ_3} & 0 \\0 & 0 & 1\end{pmatrix}`$,

   ただし、$`A`$は格子基底 $`a_1`$, $`a_2`$, $`a_3`$ を各行とする3×3行列、$`L`$ は $`LL^T = AA^T`$ を満たす下三角行列。

3. PCシフト $`Δx`$, $`Δy`$, $`Δz`$ とその推定誤差、入力ファイルdata.txt作成時に仮定したカメラ長, パターンセンターの座標を $`L^{old}`$, $`(X^{old}, Y^{old})`$ としたとき、新しいカメラ長 $`L^{new}`$ ととパターンセンターの座標 $`(X^{new}, Y^{new})`$ は以下に等しい:

   $`L^{new}=(1-Δz)L^{old}`$
   
   $`(X^{new}, Y^{new}) = (X^{old}, Y^{old}) + (L^{old}Δx, L^{old}Δy)`$

これら out.txt 内のパラメータに関する注意点として，
- まれに，$`M^{new}`$ やカイ二乗値などの指標 が精密化によって悪くなるが，これは、各バンドに割り当てられるミラー指数が精密化の結果，変わるため。
- 上記の推定誤差は， 入力角 $`φ`$, $`σ`$, $`σ_{begin}`$, $`σ_{end}`$ が 1度程度の誤差を有するという仮定の下で 非線形最小二乗法を実施したときの波及誤差になる．

### Figure_of_merit_Mnew
$`M^{new}は粉末回折のde Wolff Mの一般化として定義されるため， Mは対称性の高い格子を好むという点を除けば， 非常によく似た性質を持ちます ([詳細](https://github.com/rtomiyasu/ProjectEBSDConograph/blob/main/figures/table5_2_jp.png))．すなわち，
- ある格子定数について$`M^{new} > 10`$ なら，指数付けに成功した可能性が高い。
- 正しい格子定数は，最も大きな $`M^{new}`$ を得た格子定数のうちのどれかと考えられる。
- ほぼ同じ格子定数はBravais typeが異なる場合も，$`M^{new}`$ の値はほぼ同じになる。

したがって，リストから正しい格子定数を選ぶ際は， figure of meritの値と，どのブラベー格子に分類されているかの両方をチェックしてください．

## うまくいかないときに変更するinput_txtのパラメータ
探索方法は以下の2通りから選択できます。
- 高速探索 (対称性が高いユニットセルには十分な場合が多い)
- 徹底探索 (全てのケースに有効)

以下の値を大きくすると探索領域はさらに拡大できますが(上のパラメータほど影響力が強い, input.txt内の説明も参照)， すでに大きめの値がデフォルト値として設定されています:

- Upper bound on errors in phi, sigma, sigma_begin, sigma_end (入力角の誤差の上界): 1 度,
- Max |h|,|k|,|l| used for indexing (指数付けに用いる|h|,|k|,|l|の最大値): 7,
- Tolerance level for errors in the unit-cell scales (格子定数スケールの誤差の許容レベル): 3.,
- Resolution for Bravais-type determination (ブラベー格子決定で仮定する格子定数の解像度): 0.02.

## バグ報告のお願い
コード改善のために，疑わしい事例を見つけられましたら，以下のアドレスまでご連絡頂けると助かります．
その際，"LOG_CONOGRAPH.txt"を含む全ての入出力ファイルをご添付ください．

- tomiyasu (at) imi.kyushu-u.ac.jp

## 論文の引用に関して
EBSD-CONOGRAPHを用いて得られた研究成果を発表される際は，以下の文献を引用してください．

 R. Oishi-Tomiyasu, T. Tanaka, J. Nakagawa, “Distribution rules of systematic absence and generalized deWolff figure of merit applied to EBSD ab-initio indexing”, [web page](https://onlinelibrary.wiley.com/iucr/doi/10.1107/S1600576721002120).
