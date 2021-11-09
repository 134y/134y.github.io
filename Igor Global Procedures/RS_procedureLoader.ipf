#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma version = 2.0.0
#pragma ModuleName = RSprocLoader

#include "RS constant"
#include "RS hotfix"
#include "RS takadaShimadaExtension"
#include "RS programingTools"
#include "RS uIexpansion"
#include "RS fileLoader"
#include "RS fitFunctions"
#include "RS multiFuncCurveFit"
#include "RS procBrowser"

strconstant ksFuncParamDeclList = "variable|variable|string|wave|dfref|funcref|struct|int|int64|uint64|double|complex"

Menu "Macros"
	"RS Procedure Browser /P",/Q, ProcedureBrowser()
	"Reload Included Procedures /R",/Q, reloadProc()
	"-"
	Submenu "my procs"
		"Open procedureloader Porc;-;-",/Q, DisplayProcedure/W=$"RS_procedureLoader.ipf"
		"Load Graph Display extension", /Q, Procloader("\"RS GraphDisplay\"")
		"Load Wave Transform", /Q, Procloader("\"RS WaveTransform\"")
		"Load Wave Arithmetic", /Q, Procloader("\"RS WaveArithmetic\"")
		"Load SVD analysis", /Q, Procloader("\"RS Svd\"")
		"Load Cosmic rays remover", /Q, Procloader("\"RS CosmicRayRemove\"")
		"Load Anscombe transform", /Q, Procloader("\"RS Anscombe\"")
		"Load Infojmation analysis", /Q, Procloader("\"RS Info\"")
		"-"
		"\\M0: /L<B:Load above all at once",/Q,Procloader("\"RS GraphDisplay\"");Procloader("\"RS WaveTransform\"");Procloader("\"RS WaveArithmetic\"");Procloader("\"RS Svd\"");Procloader("\"RS CosmicRayRemove\"");Procloader("\"RS Anscombe\"");Procloader("\"RS Info\"")
		"-"
		"Load Distortion Calibration package",/Q,Procloader("\"RS DistCalibration\"")
		SubMenu "Wavelength Calibration"
			"Load" ,/Q, Procloader("\"RS WlCalibration\"")
		end
		"-"
		"Load Scaling Adjuster for Images",/Q, Procloader("\"RS ImageScale\"")
		"Load Image Display extension",/Q, Procloader("\"RS ImageDisplay\"")
		"Load AutoSize Images2",/Q, Procloader("\"Autosize Images2\"")
		"-"
		"Load Moving window Fitting",/Q, Procloader("\"RS MovingWindowFitting\"")
		"Load Baseline Estimation",/Q, Procloader("\"RS BaseLineEstimation\"")
		"Load WaveFit",/Q, Procloader("\"RS WaveFit\"")
		"Load GlobalFit",/Q, Procloader("\"RS GlobalFit\"")
		"Load SimpleHamand",/Q, Procloader("\"RS SimpleHamand\"")
		"-"
		"Load Matrix Factorization display",/Q, Procloader("\"RS MFDisplay\"")
		"Load genetic NMF analysis",/Q, Procloader("\"RS GenNMF\"")
		"Load Cluster analysis",/Q, Procloader("\"RS Cluster\"")
		"Load Discriminant analysis",/Q, Procloader("\"RS LDA\"")
		"Load PLS analysis",/Q, Procloader("\"RS PLS\"")
		"Load Hadamard analysis",/Q, Procloader("\"RS Hadamard\"")
		"-"
		"Load Spectrum Simulator",/Q, Procloader("\"RS SimulateSpectrum\"")
		"-"
		"Load Raster Image construction",/Q, Procloader("\"RS RasterImage\"")
		"Load Spiral Image construction",/Q, Procloader("\"RS SpiralImage\"")
		"Load SERDS analyzer",/Q, Procloader("\"RS SERDSanalysis\"")
		"Load SIRMS analyzer",/Q, Procloader("\"RS SIRMSAnaly\"")
		"Load SIRM analyzer",/Q, Procloader("\"RS SIRM\"")
		"Load Particle analyzer",/Q, Procloader("\"RS Particle\"")
		"Load Thermo analysis",/Q, Procloader("\"RS ThermoAnalysis\"")
		"Load NIR analyzer",/Q, Procloader("\"RS NIRanalysis\"")
		"Load fsTRIR analyzer",/Q, Procloader("\"RS fsTRIRanalysis\"")
		"Load Camera support", /Q, Procloader("\"RS CameraSupport\"")
	end
end

Menu "Macros"
	Submenu "WM procs"
		"Load matrix from table", /Q, Procloader("<MatrixFromTable>")//MatrixFromTableLoader()
		"Load IgorThief", /Q,  Procloader("<IgorThief>")//IgorThiefLoader()
		"Load CmplxToMagPhase", /Q,  Procloader("<CmplxToMagPhase>")//CmplxToMagPhaseLoader()
		"Load Decimation", /Q,  Procloader("<Decimation>")//DecimationLoader()
		"Load Manual Peak Adjust", /Q,  Procloader("<Manual Peak Adjust>")//ManualPeakAdjustLoader()
		"Load Peak AutoFind", /Q,  Procloader("<Peak AutoFind>")//PeakAutoFindLoader()
		"Load Append Calibrator", /Q,  Procloader("<Append Calibrator>")//AppendCalibratorLoader()
		"Load Straighten Tags", /Q,  Procloader("<Straighten Tags>")//StraightenTagsLoader()
		end
end

Menu "Macros"
	SubMenu "-"
	end
	subMenu "Quick"
		MenuFunctionList(), ExecuteProcMenuSelection()
	end
end

function ProcLoader(str)
	string str
	string cmd
	sprintf cmd,"INSERTINCLUDE "+str
	Execute /P/Q/Z cmd
	Execute /P/Q/Z "COMPILEPROCEDURES "
end

Function reloadProc()
	string procwinlist = winlist("RS *",";","WIN:128")
	variable num = itemsinlist(procwinlist)
	variable i
	string text
	string cmd

	if(NumberByKey("IGORVERS", IgorInfo(0))  == 7 )	// if IP7.00 , do nothing. (possibly a bug?)
		return 0
	endif
	for(i=0;i<num;i+=1)
//		sprintf cmd, "CloseProc /Name=\"%s\"/COMP",stringfromlist(i,procwinlist)
		sprintf cmd, "CloseProc /Name=\"%s\"%s",stringfromlist(i,procwinlist),selectstring(i==num-1,"","/COMP")
		execute/P/Q cmd
	endfor
	print "Procs reloaded", time()
end

function AfterCompiledHook()
	ClearCaller()
	DefaultGUIControls os9
end

Function HandleMenuTemplate(str)
	string str
	str =  stringfromlist(0,str,"(")
	if (GetkeyState(0) & 2^2)  // shift key pressed
		DisplayProcedure str
	else
		ToCommandLine ProcsTemplate(str,"")
	endif
end

Static Function/S ParseProcsText(procName, filename, [multiline])
	// Returns "Function/... funcName(.....) : subtype", or "" if the proc doesn't exist
	String procName,fileName
	variable multiline		// BOOL (0): if set, output multiline procstext as it is
	Variable lines = 0
	variable i
	string regExprStr = "(.*?(?P<pn>\(((?>[^()]+)|(?P>pn))*\)))"  // capturing outmost pairs of nested parenthesis (see Igor Help > Programming Techniques > Regular Expressions > Recursive patterns)
	string procStr, LineStr, tempStr
	SplitString/E=(regExprStr) ProcedureText(procName, lines, fileName), procStr

	if(multiline)
		return procStr
	endIf

	regExprStr = "(.*?)(//|$)"  // capture#1 line contents without the comment section
	for(i = 0; i < itemsinList(procStr,"\r"); i += 1)
		LineStr = StringFromList(i, procStr, "\r")
		SplitString/E=(regExprStr) LineStr, tempStr
		procStr = AddListItem(tempStr, RemoveListItem(i, procStr, "\r"), "\r", i)
	endfor
	return ReplaceString("\r", procStr, " ", 0, Inf)
end

Function/S ProcsTemplate(procName, fileName)
	String procName, fileName

	string moduleName = ""
	String template= ParseProcsText(procName,fileName)	// "Function/C funcName(.....) : subtype", or "" if the proc doesn't exist
	if(strlen(template))
		if (strsearch(procName, "#", 0) >= 0)  // procName contains moduleName
			moduleName = StringFromList(0, procName, "#") + "#"
			procName = StringFromList(1, procName, "#")
		endif
		Variable argsStart= strsearch(template, procName, strsearch(template, "(", 0), 3)
		Variable argsEnd= strsearch(template,")",argsStart)
		string args = template[argsStart,argsEnd]
		LogDebug({"template:", template, "args: ", args}, getCaller(), kDebug)
		string regExprStr = "(?i)^(.*?[\(|\s]+)(?:\s?(?:"+ ksFuncParamDeclList + ")(?:\s*?/\S+)*?)(\s+\w.*$)"
		string preStr = args, postStr = ""
		do	// removes parameter declaration for inline parameter syntax
			args = preStr + postStr
			SplitString/E=(regExprStr) args, preStr, postStr
		while(V_flag)
		template = moduleName + args
	endif
	template = ReplaceString("\t", template, "", 0, Inf)
	template = ReplaceString(" ", template, "", 0, Inf)
	template = ReplaceString(",", template, ", ", 0, Inf)
	return template
End

Function/S ProcsFirstLine(procname,filename)
	//* Return any string enclosed inside double quotation mark at the first line of function after the function name line.
	String procName, fileName

	Variable lines=0
	String text=ProcedureText(procName,lines, fileName)
	text= StringFromList(1,text,"\r")
	text = stringfromlist(1,text,"\"")
	if(strlen(text)==0)
		text=""
	endif
	return text
End

Function/S MenuFunctionList()
	string list =""
	string funclist = functionlist("*menu",";","Kind:2")
	variable num = itemsinlist(funclist)
	string str
	variable i

	for(i=0;i<num;i+=1)
		str = ProcsFirstLine(stringfromlist(i,funclist),"")
		if(strlen(str) == 0)
			str = "[" + stringfromlist(i,funclist) + "]"
		endif
		list = addlistitem(str,list,";",Inf)
	endfor

	return list
end

Function/S AutoMenuFuncList(procedureWinTitleStr)
	string procedureWinTitleStr
	string list =""
	string funclist = functionlist("*",";","Kind:18,WIN:"+procedureWinTitleStr)		// kind: list normal and override user-defined functions including static
	variable num = itemsinlist(funclist)
	variable i
	string infostr
	string tempstr
	string prestrList = "\\M0::;\\M0:(:  "
	variable FuncSpec

	for(i=0;i<num;i+=1)
		tempstr = ProcsFirstLine(stringfromlist(i,funclist),procedureWinTitleStr)
		if(strlen(tempstr))
			if(stringmatch(tempstr,"-*"))
				list = addlistitem(tempstr,list,";",Inf)
			else
				infostr = functioninfo(stringfromlist(i,funclist),procedureWinTitleStr)
				FuncSpec = stringmatch("static", stringbykey("SPECIAL",infostr))
				list = addlistitem(stringfromlist(FuncSpec,prestrList) + ProcsTemplate(stringfromlist(i,funclist), procedureWinTitleStr) + "	: "+tempstr,list,";",Inf)
			endif
		else
			infostr = functioninfo(stringfromlist(i,funclist),procedureWinTitleStr)
			FuncSpec = stringmatch("static", stringbykey("SPECIAL",infostr))
			list = addlistitem(stringfromlist(FuncSpec,prestrList) + ProcsTemplate(stringfromlist(i,funclist), procedureWinTitleStr),list,";",Inf)
		endif
	endfor
	return list
end

Function  ExecuteProcMenuSelection()
	GetLastUserMenuInfo
//	print S_Value, V_value, V_flag
	Execute/P/Z stringfromlist(V_value-1,functionlist("*menu",";","Kind:2"))+"()"
end

Function  ExecuteAutoMenuFuncSelection(procedureWinTitleStr)
	string procedureWinTitleStr
	GetLastUserMenuInfo
//	print S_Value, V_value, V_flag
	string fname = stringfromlist(0,S_value,"(")
	HandleMenuTemplate(fname)
end
