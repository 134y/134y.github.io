---
title: 出力先指定方法
date: 2021_0217
tags: 
---

## 概要

- 出力先を dest (wave reference) にて指定(省略可)
- dest 内容をもとに下記のいずれかのwave formで出力
	+ dest が示す既存waveに上書き
		* dest = global wave, dest = self() の時のみ
	+ 自動で新しく名前をつけ、global wave を新たに作成(既にあるなら上書き)
	+ free wave
+ Workers
	- PrepareDest: dest を解釈, 名前の自動生成など
	- StdOutput: 出力結果をdestで示す内容に実体化

```IgorPro
fuction AnalysisFunc(wave srcWav, [wave/Z dest])
	
	wave destWav = Analysis_Core_func(srcWav)
	StdOutput(destWav, PrepareDest(dest, srcWav, prefix, suffix))
end
```

### Acceptable input types of `dest`

- $""  : null input (省略時の設定でもある)
- waveref (srcWav ならoverwrite)
- name("destwavname")
- self() or name("_self_")
- name("") or free()  : 必ず free waveとして出力される

### Output types and avilable options

#### 設計

- overwrite は 基本srcWavに対して上書き(srcWavのあるfunctionは基本 srcWavを入れる)
- 同一Functionから複数wave出力がある場合の区別のため fix/lastFixを用意
	- fix/lastFix がある場合のoverwriteは上書きではなくsrcにfix/lastFixを加えたものを出力
- Default (dest = null) はかなり柔軟(場合わけ多し)な出力設定に
	- constantで挙動を変えられるようにしても良いかも...

#### `PrepareDest(dest, srcWav, pre, suf)`

| dest\scrWav |      global       | name("")/free  |     name(str)     |       null        |
|:-----------:|:-----------------:|:--------------:|:-----------------:|:-----------------:|
|   usecase   |       user        |      user      | preset(no srcWav) |         ?         |
|    -----    |       -----       |     -----      |       -----       |       -----       |
|    null     | name(pre+src+suf) |  **name("")**  | name(pre+str+suf) | **name(pre+suf)** |
|  name(fmt)  |  name(fmt;@=src)  | name(fmt;@="") |  name(fmt;@=str)  |  name(fmt;@="")   |
|   self()    |      srcWav       |     srcWav     |      srcWav       |     name("")      |
|   free()    |     name("")      |    name("")    |     name("")      |     name("")      |
|    free     |     name("")      |    name("")    |     name("")      |     name("")      |
|   global    |       dest        |      dest      |       dest        |       dest        |

src = NameOfWave(srcWav)
srcWav: srcWavに入れたreferenceそのもの (以前のOverwrite)
dest: destに入れたglobal referenceそのもの (以前のOverwrite)

##### 挙動

- name(fmt) で決めうち文字列の場合、必ずその文字列で出力される
	+ **同じFunction内でこのstyleの`PrepareDest`が複数あるとuser inputによっては結果が重ね書きされる可能性**
- User input によって **srcWavのOverwrite が可能**
+ srcWav = free/null, dest = null の場合 (入力が前処理の解析結果流し込みなど)
	+ 初版: name(pre+suf)
		+ pre/suf がemptyかどうかで global or free どちらの出力にもなりうる。これは問題...
		+ 一律 free にしてしまうのが運用上楽? (正しい名前の推定は難しい)
	+ 改訂版: **src = empty ならfree wave出力, null は(fix使用時との整合性のため) global出力することに**
- srcWav = null は特殊な状況, user の inputとしてはなさそう、fix を使うときとの整合性のため...との位置づけ
- dest/srcWav共に free()と一般のfree waveは全く同じ扱い

##### 用途

- 1つしかないorメインの出力結果

#### `PrepareDest(dest, srcWav, "", "", fix = "fix", lastFix = "lf")`

| dest\scrWav |         global         |     name("")/free     |       name(str)        |         null          |
|:-----------:|:----------------------:|:---------------------:|:----------------------:|:---------------------:|
|   usecase   |          user          |         user          |           ?            |   preset(no srcWav)   |
|    -----    |         -----          |         -----         |         -----          |         -----         |
|    null     |    name(fix+src+lf)    |       name("")        |    name(fix+str+lf)    |     name(fix+lf)      |
|  name(fmt)  | name(fix+fmt+lf;@=src) | name(fix+fmt+lf;@="") | name(fix+fmt+lf;@=str) | name(fix+fmt+lf;@="") |
|   self()    |    name(fix+src+lf)    |     name(fix+lf)      |    name(fix+str+lf)    |     name(fix+lf)      |
|   free()    |        name("")        |       name("")        |        name("")        |       name("")        |
|    free     |        name("")        |       name("")        |        name("")        |       name("")        |
|   global    |    name(fix+de+lf)     |    name(fix+de+lf)    |    name(fix+de+lf)     |    name(fix+de+lf)    |

src = NameOfWave(srcWav)
de = NameOfWave(de)

##### 挙動

- name(fmt) で決めうち文字列の場合でもプログラム的に必ず文字列"fix/lf"を付加できる
	+ **同じFunction内で複数の出力結果を作成する場合にこれを利用**
- srcWavのOverwriteを防げる
+ srcWav = free/null, dest = null の場合
		+ 初版: fix+lf で global waveが 必ず作られる
		+ 一律freeでも良い?? だめ 
		+ 改訂版: **src = empty ならfree wave出力, null は(fix使用時との生合成のため) global出力することに**
- dest free()と一般freeは全く同じ扱い

##### 用途

- メイン出力と出力先が重複しないことを優先させたいsub出力用
- srcWavと異なるwave構造を持つwaves

#### Examples

標準input : srcWav = global, dest = null とする

##### srcWave名base

+ `PrepareDest(dest, srcWav, prefix, suffix)` vs `PrepareDest(dest, srcWav, "", "", lastFix = "_suffix"))`
	- 標準input	=>	name(src_suffix)
	- srcWav = free/free() , dest= null	=> どちらもfree()出力
	- name(srt)	=>	どちらも @ = src
	- self()	=>	前者はsrcWav 上書き, 後者は default出力
	
##### 固定名baseの使い分けについて

+ `PrepareDest(dest, srcWav, "M_GenAnsFCoef", "", format = "<"))` vs `PrepareDest(dest, name("M_GenAnsFCoef"), "", ""))`
	- 標準input	=>	どちらも name("M_GenAnsFCoef") 
	- srcWav = free/free() , dest= null	=>	前者はfree()出力 , 後者はglobal出力
	- name(srt)	=>	前者は @ = src, 後者は @ = "M_GenAnsFCoef"
	- self()	=>	前者はsrcWav 上書き, 後者は default出力
+ `PrepareDest(dest, srcWav, "", "", fix = "M_baseWav", format = "")` vs `PrepareDest(dest, $"", "", "", fix = "M_baseWav")`
	- 標準input => どちらも name("M_baseWav") 
	- srcWav = free/free() , dest= null	=>	前者はfree()出力 , 後者はglobal出力
	- name(srt)	=>	前者は @ = src, 後者は @ = ""
	- self()	=>	どちらも name("M_baseWav") 

よって...

- srcWave あり、 1つ or メイン出力: 
	- 元のwave名base: `PrepareDest(dest, srcWav, "", "_suffix"))`
	- 固定名base: `PrepareDest(dest, srcWav, "W_newname", "", format = "<>"))`
- srcWave なし、 1つ or メイン出力: 
	- 元のwave名base: N/A
	- 固定名base: `PrepareDest(dest, name("W_newname"), "", ""))`
- srcWave あり、 複数 or サブ出力: 
	- 元のwave名base: `PrepareDest(dest, srcWav, "", "", lastFix = "_suffix"))`
	- 固定名base: `PrepareDest(dest, srcWav, "", "", fix = "M_baseWav", format = "<>")`
- srcWave なし、 複数 or サブ出力: 
	- 元のwave名base: N/A
	- 固定名base: `PrepareDest(dest, $"", "", "", fix = "M_baseWav")`




free + "" なら freewav
free + "destwavname" なら duplicate/O destWav, $destwavname
free + "_self_" なら duplicate/O destWav, srcWav
normal なら duplicate/O destWav, dest
null なら duplicate/O destWav, $defaultname

- (default) normal wave 生成 (名前も自動でつける)
	- dest = $"" 
	- (存在しない実在waveを誤って指定した場合もこれ)
		- WaveType(dest, 2) == 0 // null
		- destWavname = "<@>"
		- wave retWav = name("<@>")
- normal wave 生成 (名前指定)
	- dest = name("destwavname")
		- WaveType(dest, 2) == 2 // free wave
		- strlen(destwavname) > 0
		- wave retWav = dest = name("destWavname")
- 実在のwave(に上書き)
	- dest = waveref (srcWav ならoverwrite)
		- wave retWav = dest
- srcWaveに上書き (overwrite)
	- dest = self() or $self
		- WaveType(dest, 2) == 2 // free wave
		- CmpStr(destwavname, "_self_") == 0  // equal
		- wave retWav = dest
- free wave (名付け不要)
	- dest = name("") or free() or $free
		- WaveType(dest, 2) == 2 // free wave
		- strlen(destwavname) == 0


StdOutput処理

dest, destWav, destWavname


strlen(destWavname) == 0

duplicate/O destWav, retWav



(最後は duplicate/O のターゲットを destにする)
上書きはなるだけ最後 (error終了時に元のwaveを保護するため)



real wave 
-> accidental overwrite
-> with name
-> no-auto deletion

free wave 
-> no accidental overwrite
-> no name
-> chance of memory leak
-> auto deletion

