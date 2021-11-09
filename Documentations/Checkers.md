---
title: Argument checkについて  
date: 
tags: 
---


CheckWave(condition, enforce)

Action
	condition (条件に適合しなければerror)
	enforce (条件に適合しなければ強制させる) NAN処理とdimくらいか


Wave property

exsits: BOOL: WaveExisits()
type: "null; global; free" : WaveType(, 2)
dim: 0-4: WaveDims()
dims: {,,,}: DimSize()
numTypeNormal: BOOL: numType(sum()) normal; NAN/±inf

Waves compare

