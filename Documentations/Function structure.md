
Function WaveAnalysis(													// Initialization
	wave srcWav,
	variable arg1,
	string arg2,
	[wave/Z dest, skipWaveCheck])

	STRUCT FuncBreadcrumb fb
	[fb] = openFuncBreadcrumb("SVD", {srcWav}, ParamIsDefault(verbose) ? kVerbose : verbose)
	ArgCheck()											// Validate input waves/args
	// Initialize variables
	// 

	core logic


	//RegenCommand
	wave destMain = BuildDest(dest, srcWav, "", "prefix_", "_suffix")	// Prepare destination waveref
	wave retWav = StdOutput(destWav, destMain)							// generate destination wave
	output

end

Structure FuncTracking
	variable timer
	string funcname 
	string caller
endstructure

mulitple Return syntaxはコマンドラインから使えない => Base functionでは用いないことにする


### dest導入に際しての一括変換

(summarize\(.*?)\s?destwavname\s?=\s?\"\"\s?(.*?\))
$1 dest = free()$2



### FuncBreadcrumb導入に際しての一括変換

search
	variable timer = ticks
	string funcname = \"(\S*?)\"
	string caller = setCaller\(\"\[\"\s?\+\s?funcname\s?\+\s?\"\]\"\)
	verbose = selectnumber\(paramisdefault\(verbose\),\s?verbose,\s?kVerbose\)
replace
	STRUCT FuncBreadcrumb fb
	[fb] = openFuncBreadcrumb("$1", {srcWav}, ParamIsDefault(verbose) ? kVerbose : verbose)

search
(//)?	LogInfo\(\{\"Elapsed time: ", Secs2Time\(\(ticks-timer\)/60.15,5\)\}, caller, verbose\)([\s\S]*?)
	setCaller\(\"\"\)
replace
$1$2	closeFuncBreadcrumb(fb)

search
(log\S*?\(.*?(?<!fb)\W)(verbose|caller)(.*?\))
replace
$1fb.$2$3

search -> replace
destWav, funcname -> destWav, fb.funcname
\$funcname -> $fb.funcname

fb.に続かないverboseの検索
(?<!fb)\W(caller|verbose)
