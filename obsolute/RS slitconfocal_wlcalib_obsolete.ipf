#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma version = 2.0.0
#pragma ModuleName = RSSCwlCalib

#include "RS GlobalFit"
#include "RS WLCalibration"
#include "RS WaveTransform"
#include "RS WaveArithmetic"
#include <Peak AutoFind>

Menu "Macros"
	SubMenu "-RS"
	end
	SubMenu "RS slitconfocal_wlcalib"
		"<Ishow procedrure;-",/Q, DisplayProcedure/W=$"RS slitconfocal_wlcalib.ipf"
		AutoMenuFuncList("RS slitconfocal_wlcalib.ipf"),/Q,ExecuteAutoMenuFuncSelection("RS slitconfocal_wlcalib.ipf")
	end
end

function QMseparatorRSSCwlCalibMenu()
//		"-"		//		"-;(<IRSSCwlCalib;-"
end

Menu "Calibration"
	 	"\\M0::"+ProcsTemplate("Autofit",""),/Q, HandleMenuTemplate("autofit(refwav,[peak_obs,peak_ref,ROI,polyorder])")
 		"\\M0::"+ProcsTemplate("MakeNewWlq",""),/Q, HandleMenuTemplate("MakeNewWlq(paramwav,qnum,[xdim])")
 		"\\M0::"+ProcsTemplate("AlignAbscissa",""),/Q, HandleMenuTemplate("AlignAbscissa(refwav,paramwav,xdestwav,[mode,exwl])")
 		"\\M0::"+ProcsTemplate("generateXwav",""),/Q, HandleMenuTemplate("generateXwav(newXwavname,startX,endX,inc)")
 		"\\M0::"+ProcsTemplate("AlignBothAxes",""),/Q, HandleMenuTemplate("AlignBothAxes(refwav,xparamwav,xdestwav,yparamwav,ydestwav,[exwl])")
end


function autofit(refwav,[paramDFR,peak_obs,peak_ref,ROI,polyorder,exwl,q1,q2, smoothing, maskwav, centroid, suffix, verbose])
	wave refwav
	wave/Z peak_obs,peak_ref
	variable ROI
	variable polyorder
	variable exwl  // only used to calculate the RamanShift range
	variable smoothing
	variable centroid
	variable q1,q2
	wave maskwav
	DFREF paramDFR  // if any of the contents are specified by their individual option parameter, those are given priority over the waves in the paramDFR
	string suffix
	variable verbose

	variable timer = ticks
	string funcname = "autofit"
	string caller = setCaller("[" + funcname + "]")
	suffix = selectstring(paramisdefault(suffix),suffix,"")
	verbose = selectnumber(paramisdefault(verbose), verbose, kVerbose)

	string paramwavname = ValidateName("M_WlCalibParam" + suffix, overwrite = 1)
	if(cmpstr(paramwavname, ksValidateNameCanceled) == 0)
		setCaller("")
		return 0  // user canceled
	endif
	string paramsigmawavname = ValidateName("M_WlCalibParam_sigma" + suffix, overwrite = 1)
	if(cmpstr(paramsigmawavname, ksValidateNameCanceled) == 0)
		setCaller("")
		return 0  // user canceled
	endif

	wave xwav = getXwav(refwav)

	// loading the
	if(!DataFolderRefStatus(paramDFR)) // param is invalid
		DFREF paramDFR = root:Packages:RSwlcalib
	endif
	if(!WaveExists(peak_obs))
		wave peak_obs = $(getDataFolder(1,paramDFR) + "parameters:peak_obs")
	endif
	if(!WaveExists(peak_ref))
		wave peak_ref = $(getDataFolder(1,paramDFR) + "parameters:peak_ref")
	endif
	if(ParamIsDefault(ROI))
		NVAR ROIwidth = $(getDataFolder(1,paramDFR) + "parameters:ROIwidth")
		ROI = ROIwidth
	endif
	if(ParamIsDefault(polyorder))
		NVAR gpolyorder = $(getDataFolder(1,paramDFR) + "parameters:polyorder")
		polyorder = gpolyorder
	endif
	if(ParamIsDefault(exwl))
		NVAR gexwl = $(getDataFolder(1,paramDFR) + "parameters:exwl")
		exwl = gexwl
	endif
	// wave lam = $(getDataFolder(1,paramDFR) + "parameters:lam")
	// variable lam_offset = str2num(GetNoteProc(lam,"RecalcLam","offset"))
	// variable lam_slope = str2num(GetNoteProc(lam,"RecalcLam","slope"))
	// variable lam_xdim = str2num(GetNoteProc(lam,"RecalcLam","xdim"))
	// duplicate/O peak_ref peak_ref_pixel
	// peak_ref_pixel = (peak_ref[p] - lam_offset) / lam_slope + lam_xdim / 2
	wave peak_ref = peak_obs
	q1 = selectnumber(ParamIsDefault(q1),q1,0)
	q2 = selectnumber(ParamIsDefault(q2),q2,dimsize(refwav,1) - 1)

	LogInfo({"peak_obs: ", varwave2list(peak_obs,listsepstr=",")}, caller, verbose)
	LogInfo({"peak_ref: ", varwave2list(peak_ref,listsepstr=",")}, caller, verbose)
	LogInfo({"ROI: ", num2str(ROI)}, caller, verbose)
	LogInfo({"polyorder: ", num2str(polyorder)}, caller, verbose)
	LogInfo({"exwl: ", num2str(exwl)}, caller, verbose)
	LogInfo({"q1: ", num2str(q1)},caller,verbose)
	LogInfo({"q2: ", num2str(q2)},caller,verbose)

//	print peak_obs
	variable xdim = dimsize(refwav,0)
	variable ydim = dimsize(refwav,1)
	make/n=4/FREE/D dims = dimsize(refwav, p)
	LogInfo({"dims:",varwave2list(dims,listsepstr=",")}, caller, verbose)
	variable numpeak = dimsize(peak_obs,0)

	variable i, j, l, m
	variable pcsrA,pcsrB
	variable V_fitoptions = 4
	variable V_FitError = 0
	string win

	make/n=(dims[1],polyorder+1)/O $paramwavname/wave=M_WlCalibParam
	make/n=(dims[1],polyorder+1)/O $paramsigmawavname/wave=M_WlCalibParam_sigma
	make/n=0/O W_chisq
//	edit M_WlCalibParam

	make/n=(dims[0])/O temp
//	duplicate/FREE peak_obs
	make/n=(dims[1],numpeak,dims[2])/O M_temp_peakposi,M_temp_peakposi_sigma,fit_M_temp_peakposi
	make/n=(dims[1])/FREE startwl, endwl
	Make/O/T/N=1/FREE T_Constraints
	T_Constraints[0] = {"K1 > 0","K2>dummy","K2<dummy","K3<"+num2str(2*ROI),"K3>0"}
	M_temp_peakposi = NAN

	m = 0
	do
		for(i=0;i<dims[1];i+=1)

			temp = refwav[p][i][m]

			for(j=0;j<numpeak;j+=1)

				if(numtype(peak_obs[j]) != 0)
					continue
				endif

				if(centroid == 0)
					pcsrA = peak_obs[j] - ROI
					pcsrB = peak_obs[j] + ROI
					T_Constraints[1] = "K2>"+num2str(pcsrA)
					T_Constraints[2] = "K2<"+num2str(pcsrB)

					V_FitError = 0
					try
						curvefit/Q/W=2/N=1 gauss, temp[pcsrA,pcsrB] /C=T_Constraints//; AbortOnRTE
					catch
						if (V_AbortCode == -4)
							Variable CFerror = GetRTError(1)	// 1 to clear the error
							LogError({"curve fit: ", GetErrMessage(CFerror)}, caller)
							LogDebug({"i: ", num2str(i),"\tj: ", num2str(j)}, caller, kDebug)
						endif
					endtry
					wave w_sigma,w_coef

					if(V_FitError || k2 <pcsrA || K2 > pcsrB || w_sigma[2] > 1)		// Fitting Error
						M_temp_peakposi[i][j][m] = NAN
						M_temp_peakposi_sigma[i][j][m] = NAN
						V_FitError = 0
					else
						M_temp_peakposi[i][j][m] = w_coef[2]
						M_temp_peakposi_sigma[i][j][m] = W_sigma[2]
					endif
				else
					for(l = ROI; l < ROI * 4; l += 1)
						if(wavemin(temp,peak_obs[j] - l, peak_obs[j]) == wavemin(temp, peak_obs[j] - l - 1, peak_obs[j]))
							break
						endif
					endfor
					pcsrA = peak_obs[j] - l
					for(l = ROI; l < ROI * 4; l += 1)
						if(wavemin(temp, peak_obs[j], peak_obs[j] + l) == wavemin(temp, peak_obs[j], peak_obs[j] + l + 1))
							break
						endif
					endfor
					pcsrB = peak_obs[j] + l

					// display/N=tempgraph/K=1 temp
					// ModifyGraph mode=4, marker=19
					// setaxis bottom, pcsrA - 2, pcsrB + 2
					// SetAxis/A=2 left
					// Tag/C/N=text0 temp, pcsrA,"\\OX"
					// Tag/C/N=text1 temp, pcsrB,"\\OX"
					// PauseForUser tempgraph

					bandStatsXY(getXwav(temp),temp,pcsrA,pcsrB)
					NVAR V_BandCentroid
					M_temp_peakposi[i][j] = V_BandCentroid
					bandStatsXY(getXwav(temp),temp,pcsrA - 1,pcsrB + 1)
					NVAR V_BandCentroid
					M_temp_peakposi_sigma[i][j] = abs(M_temp_peakposi[i][j] - V_BandCentroid)/sqrt(2)
				endif
			endfor	// j
		endfor	// i
		m += 1
	while(m < dims[2])

	if(smoothing)
		if(!WaveExists(maskwav))
			make/N=(dims[1])/FREE maskwav = NAN
			maskwav[q1,q2] = 1
		endif
		for(j=0;j<numpeak;j+=1)
			if(smoothing == 1)
				make/O/N=(numpeak,smoothing + 1) $("fit_M_temp_peakposi_coef" + suffix)/wave=fit_M_temp_peakposi_coef
				CurveFit/Q/ODR=2/N=1 line, M_temp_peakposi[][j] /I=1/W=M_temp_peakposi_sigma[][j]/M=maskwav
//				CurveFit/Q/ODR=2/N=1 line, M_temp_peakposi[q1,q2][j]
				wave w_coef
				fit_M_temp_peakposi[][j] = w_coef[0] + w_coef[1] * p
			elseif(smoothing > 1)
				make/O/N=(numpeak,smoothing + 1) $("fit_M_temp_peakposi_coef" + suffix)/wave=fit_M_temp_peakposi_coef
				CurveFit/Q/W=0/N=1 poly (smoothing+1),M_temp_peakposi[][j] /I=1/W=M_temp_peakposi_sigma[][j]/M=maskwav
//				CurveFit/Q/W=0/N=1 poly (smoothing+1),M_temp_peakposi[q1,q2][j]
				wave w_coef
				fit_M_temp_peakposi[][j] = poly(W_coef,p)
			elseif(smoothing == -1) // sine function
				curvefit/Q/ODR=2 sin, M_temp_peakposi[q1,q2][j] /I=1/W=M_temp_peakposi_sigma[q1,q2][j]
				wave w_coef

				generatelinkwav(5,numpeak,linkrowlist="1;2;3;4")
				wave linkwav
				globalfit(M_temp_peakposi, skewedsin, {100,w_coef[1],w_coef[2],w_coef[3],0}, linkwav, holdstr="00000",sigmawav= M_temp_peakposi_sigma,maskwav=maskwav)
				wave w_coef
				redimension/n=(5,numpeak) w_coef
				matrixtranspose w_coef
//				fit_M_temp_peakposi[][] = mysin({w_coef[q][0],w_coef[q][1],w_coef[q][2],w_coef[q][3]},p)
				duplicate/O w_coef $("fit_M_temp_peakposi_coef" + suffix)/wave=fit_M_temp_peakposi_coef
				break
			elseif(smoothing == -2) // error function
				curvefit/Q/ODR=2 sigmoid, M_temp_peakposi[q1,q2][j] /I=1/W=M_temp_peakposi_sigma[q1,q2][j]
				wave w_coef

				generatelinkwav(4,numpeak,linkrowlist="1;2;3")
				wave linkwav
				globalfit(M_temp_peakposi, stepCONVgaussAmp, {100,w_coef[1],w_coef[2],w_coef[3]}, linkwav, holdstr="0000",sigmawav= M_temp_peakposi_sigma,maskwav=maskwav)
				wave w_coef
				redimension/n=(4,numpeak) w_coef
				matrixtranspose w_coef
//				fit_M_temp_peakposi[][] = mysin({w_coef[q][0],w_coef[q][1],w_coef[q][2],w_coef[q][3]},p)
				duplicate/O w_coef $("fit_M_temp_peakposi_coef" + suffix)/wave=fit_M_temp_peakposi_coef
				break
			endif
			fit_M_temp_peakposi_coef[j][] = w_coef[q]
			setNoteProc(fit_M_temp_peakposi,"autofit","w_coef",varwave2list(w_coef),newline=1)
		endfor

		matrixop/O res_M_temp_peakposi =  M_temp_peakposi - fit_M_temp_peakposi
		for(j=0;j<numpeak;j+=1)
			win = "autofit_peak_comp"+num2str(j)
			dowindow/F $win
			if(V_flag)
				dowindow/K $win
			endif
			display/K=1/N=$win M_temp_peakposi[][j],fit_M_temp_peakposi[][j]
			ErrorBars M_temp_peakposi Y,wave=(M_temp_peakposi_sigma[*][j],M_temp_peakposi_sigma[*][j])
			appendtograph/W=$win/L=Res_Left res_M_temp_peakposi[][j]
			res_M_temp_peakposi *= maskwav[p]
			setAxis left, wavemin(M_temp_peakposi,j*dims[1],(j+1)*dims[1]-1),wavemax(M_temp_peakposi,j*dims[1],(j+1)*dims[1]-1)
//			print wavemin(M_temp_peakposi,j*dims[1],(j+1)*dims[1]-1),wavemax(M_temp_peakposi,j*dims[1],(j+1)*dims[1]-1)
			ModifyGraph axisEnab(Left)={0,0.75}
			ModifyGraph axisEnab(Left)={0,0.75}
			ModifyGraph axisEnab(Res_Left)={0.8,1},freePos(Res_Left)=0,zero(Res_Left)=1
			ModifyGraph rgb(fit_M_temp_peakposi)=(0,0,65535)
			setaxis bottom, q1,q2
			setAxis/A=2 Res_left
			doupdate
		endfor
		tilerecentwindows(numpeak)
	else
		if(dims[2] == 0)
			duplicate/O M_temp_peakposi, fit_M_temp_peakposi,res_M_temp_peakposi
		else
			Reshape2D(M_temp_peakposi, {0}, verbose = verbose)
			duplicate/O M_temp_peakposi, fit_M_temp_peakposi,res_M_temp_peakposi
			duplicate/O peak_ref peak_ref_temp
			redimension/n=(numpeak * dims[2]) peak_ref_temp
			peak_ref_temp = peak_ref[mod(p, numpeak)]
			wave peak_ref = peak_ref_temp
		endif
		res_M_temp_peakposi=0
	endif

	win = "pause"
	dowindow/F $win
	if(V_flag)
		dowindow/K $win
	endif
	make/N=(dimsize(fit_M_temp_peakposi,1))/O peak_obs_temp
	display/K=1/N=$win peak_ref vs peak_obs_temp
	if(smoothing)
		tilerecentwindows(1+numpeak)
	else
		tilerecentwindows(1)
	endif

	for(i=0;i<dims[1];i+=1)
		peak_obs_temp = fit_M_temp_peakposi[i][p]
		ModifyGraph mode=3,marker=41
		if(polyorder == 1)
			CurveFit/Q/ODR=2/N=1 line peak_ref /X=fit_M_temp_peakposi[i][]/D/R
		else
			CurveFit/Q/ODR=2/N=1 poly (polyorder+1), peak_ref /X=fit_M_temp_peakposi[i][]/D/R
		endif
		ModifyGraph mode(fit_peak_obs)=0
		ModifyGraph zero(Res_Left)=1
//		tilerecentwindows(1)
		wave w_coef,w_sigma
		V_fiterror = 0
		addtotail(W_chisq,V_chisq)

		M_WlCalibParam[i][] = w_coef[q]
		M_WlCalibParam_sigma[i][] = w_sigma[q]

		startwl[i] = poly(w_coef,0)
		endwl[i] = poly(w_coef,dims[0])
		doupdate
		// pauseforuser pause
	endfor	// i

	LogInfo({"startwl:", num2str(wavemax(startwl)), num2str(wavemin(startwl))}, caller, verbose)
	LogInfo({"endwl:", num2str(wavemax(endwl)), num2str(wavemin(endwl))}, caller, verbose)
	LogInfo({"startwn:", num2str(10^7/exwl - 10^7/wavemax(startwl)), num2str(10^7/exwl - 10^7/wavemin(startwl))}, caller, verbose)
	LogInfo({"endwn:", num2str(10^7/exwl - 10^7/wavemax(endwl)), num2str(10^7/exwl - 10^7/wavemin(endwl))}, caller, verbose)
	LogInfo({"pixel separation (wl):", num2str((wavemax(startwl) - wavemin(endwl))/dims[0])}, caller, verbose)
	LogInfo({"pixel separation (wn):", num2str(((10^7/exwl - 10^7/wavemax(startwl)) - (10^7/exwl - 10^7/wavemin(endwl)))/dims[0])}, caller, verbose)

	note/K M_WlCalibParam, note(refwav)
		setNoteProc(M_WlCalibParam,"autofit","refwav",getWavesDataFolder(refwav,2))
		setNoteProc(M_WlCalibParam,"autofit","peak_obs",getWavesDataFolder(peak_obs,2))
		setNoteProc(M_WlCalibParam,"autofit","peak_ref",getWavesDataFolder(peak_ref,2))
		setNoteProc(M_WlCalibParam,"autofit","ROI",num2str(ROI))
		setNoteProc(M_WlCalibParam,"autofit","polyorder",num2str(polyorder))
		setNoteProc(M_WlCalibParam,"autofit","smoothing",num2str(smoothing))
		setNoteProc(M_WlCalibParam,"autofit","q1",num2str(q1))
		setNoteProc(M_WlCalibParam,"autofit","q2",num2str(q2))

	note/K M_WlCalibParam_sigma, note(M_WlCalibParam)
	note/K M_temp_peakposi, note(M_WlCalibParam)
	note/K M_temp_peakposi_sigma, note(M_WlCalibParam)
//	note/K fit_M_temp_peakposi, note(M_WlCalibParam)
	LogInfo({"Elapsed time: ", Secs2Time((ticks-timer)/60.15,5,2)}, caller, verbose)
	setCaller("")
end

Function stepCONVgaussAmp(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = f_RevstepConvGauss(x, FWHM, y0, x0) + f_StepConvGauss(x, FWHM, y1, x0)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = y0
	//CurveFitDialog/ w[1] = Amp
	//CurveFitDialog/ w[2] = x0
	//CurveFitDialog/ w[3] = FWHM

	return w[0] + f_StepConvGauss(x, w[3], w[1], w[2])
End


function MakeNewWlQ(paramwav,qnum,[xdim])
	// paramwav 	: output wave "M_WlCalibParam" from autofit function
	// qnum 		: pixel number in ydim at which the new wl wave is to be calculated
	// xdim		: number of pixels in xdimension
	wave paramwav
	variable qnum
	variable xdim
	xdim = selectnumber(ParamIsDefault(xdim),xdim,1024)

	string wlwavname = "wl_q"+num2str(qnum)

	make/n=(xdim)/O $wlwavname/wave=wlwav

	make/n=(dimsize(paramwav,1))/FREE coefwav

	coefwav = paramwav[qnum][p]

	wlwav = poly(coefwav,p)
end


function AlignAbscissa(srcWav,paramwav,xdestwav,[exwl,exwnMod, interpmode, suffix, overwrite, verbose])
	wave srcWav
	wave paramwav		//	paramwav should alwways be calculated in wavelength
	wave xdestwav
	variable exwl		// excitation wavelength (when 0, alignment will be done in wavelength instead of Raman shift)
	wave exwnMod		// wave which contains the magunitude of excitation wave modulation in wavenumber
	variable interpmode  //t=1 for linear,	t=2 (default) for cubic spline,	t=3 for smoothing spline.
	string suffix
	variable overwrite, verbose

	exwl = selectnumber(ParamIsDefault(exwl),exwl,532.2)
	interpmode = selectnumber(paramisdefault(interpmode), interpmode, 2)
	string funcname = "AlignAbscissa"
	string caller = setCaller("[" + funcname + "]")
	verbose = selectnumber(paramisdefault(verbose), verbose, kVerbose)
	suffix = selectstring(ParamIsDefault(suffix),suffix,"_L")
	suffix = selectstring(overwrite, suffix, "")
	string destwavname = ValidateName(nameofwave(srcWav) + suffix,overwrite = 1)
	if(cmpstr(destWavname, ksValidateNameCanceled) == 0)
		setCaller("")
		return 0  // user canceled
	endif

	make/n=4/FREE/D dims = dimsize(srcWav,p)
	variable xdimdest = dimsize(xdestwav,0)

	if(!WaveExists(exwnMod))
		make/FREE/N=(dims[1]) exwnMod = 0
	endif
	duplicate/FREE srcWav destwav
	redimension/n=(xdimdest,-1,-1,-1) destwav
	make/n=(dims[0])/FREE refwav_temp,xwav,dxwav
	make/n=(xdimdest)/FREE dest_temp
	make/n=(dimsize(paramwav,1))/FREE paramwav_temp
	variable i,j,l

	if(dims[1] != dimsize(paramwav,0))
		LogError({"Dimension mismatch!"}, caller)
		LogDebug({"dims[1]: ", num2str(dims[1]), "\tdimsize(paramwav,0): ", num2str(dimsize(paramwav,0))}, caller, kDebug)
		setCaller("")
		return -1
	endif

	l=0
	do
		j=0
		do
			for(i=0;i<dims[1];i+=1)
				refwav_temp = srcWav[p][i][j][l]
				paramwav_temp = paramwav[i][p]
				if(exwl)
					xwav = 10^7/exwl + exwnMod[i] - 10^7/poly(paramwav_temp,p)		// xwav in Ramanshift
				else
					xwav = poly(paramwav_temp,p)		// xwav in wavelength
				endif
				dxwav = abs(poly(paramwav_temp,p + .5) - poly(paramwav_temp,p - .5))
				variable dxmedian = median(dxwav)
				dxwav /= dxmedian
				duplicate/FREE refwav_temp, refwav_temp_didlam
				refwav_temp_didlam = refwav_temp[p] / dxwav[p]
				Interpolate2/T=(interpmode)/I=3/E=2/Y=dest_temp/X=xdestwav xwav, refwav_temp_didlam

				destwav[][i][j][l] = dest_temp[p]

			endfor			// i
			j=j+1
		while(j<dims[2])	// j
		l=l+1
	while(l<dims[3])		// l

	note/K destwav, note(srcWav)
	setNoteHistory(destwav,funcname,"srcWav",getWavesDataFolder(srcWav,2),newline=1)
	setNoteHistory(destwav,funcname,"paramwav",getWavesDataFolder(paramwav,2))
	setNoteHistory(destwav,funcname,"xdestwav",getWavesDataFolder(xdestwav,2))
	setNoteHistory(destwav,funcname,"exwl",num2str(exwl))
	setNoteHistory(destwav,funcname,"exwnMod",getWavesDataFolder(exwnMod,2))
	setNoteHistory(destwav,funcname,"interpmode",num2str(interpmode))
	setNoteHistory(destwav,funcname,"suffix",suffix, isStr = 1)
	setNoteHistory(destwav,funcname,"overwrite",num2str(overwrite))

	if(overwrite)
		duplicate/O destwav, srcWav
		wave destwav = srcWav
	else
		Duplicate/O destwav, $destwavname
		wave destwav = $destwavname
	endif
	string/g S_wavename = getWavesDataFolder(destwav,2)
	setCaller("")
end

function/WAVE AlignBothAxes(srcWav,xparamwav,xdestwav,yparamwav,ydestwav,[exwl, exwnMod, interpmode, dest, verbose])
	wave srcWav
	wave xparamwav
	wave xdestwav
	wave yparamwav
	wave ydestwav
	variable exwl		// excitation wavelength (when 0, alignment will be done in wavelength instead of Raman shift)
	wave exwnMod		// wave which contains the magunitude of excitation wave modulation in wavenumber
	variable interpmode  //t=1 for linear,	t=2 (default) for cubic spline,	t=3 for smoothing spline.
	wave/Z dest
	variable verbose

	variable timer = ticks
	string funcname = "AlignBothAxes"
	string caller = setCaller("[" + funcname + "]")
	verbose = selectnumber(paramisdefault(verbose), verbose, kVerbose)

	exwl = selectnumber(ParamIsDefault(exwl),exwl,532.2)
	interpmode = selectnumber(paramisdefault(interpmode), interpmode, 2)
	if(!WaveExists(exwnMod))
		make/FREE/N=(dimsize(ydestwav,0)) exwnMod
	endif

	duplicate/FREE/O srcWav destWav
	VolumeTranspose(destWav, {1,0,2,3}, dest = self(), verbose = verbose)
	AlignAbscissa(destWav,yparamwav,ydestwav,exwl=0, overwrite = 1, verbose = verbose)
	VolumeTranspose(destWav, {1,0,2,3}, dest = self(), verbose = verbose)
	AlignAbscissa(destWav,xparamwav,xdestwav,exwl=exwl, exwnMod=exwnMod, interpmode = interpmode, overwrite = 1, verbose = verbose)

	note/K destWav, note(srcWav)
	setNoteHistory(destWav, funcname, "srcWav", getWavesDataFolderEx(srcWav), newline = 1)
	setNoteHistory(destWav, funcname, "xparamwav", getWavesDataFolderEx(xparamwav))
	setNoteHistory(destWav, funcname, "xdestwav", getWavesDataFolderEx(xdestwav))
	setNoteHistory(destWav, funcname, "yparamwav", getWavesDataFolderEx(yparamwav))
	setNoteHistory(destWav, funcname, "ydestwav", getWavesDataFolderEx(ydestwav))
	setNoteHistory(destwav, funcname,"exwl",num2str(exwl))
	setNoteHistory(destwav, funcname,"exwnMod",getWavesDataFolderEx(exwnMod))
	setNoteHistory(destwav, funcname,"interpmode",num2str(interpmode))
	setNoteHistory(destWav, funcname, "dest", getWavesDataFolderEx(dest))

	wave retWav = StdOutput(destWav, PrepareDest(dest, srcWav, "", "_LL"))
	LogInfo({"Elapsed time: ", Secs2Time((ticks-timer)/60.15,5)}, caller, verbose)
	setCaller("")
	return retWav
end

Function AlignOrdinate(wav,pixel_old, pixel_new)
	wave wav
	wave pixel_old
	wave pixel_new

	variable interpmode = 2
	make/n=4/FREE/D dims = dimsize(wav,p)

	string ydestwavname = nameofwave(wav) + "_LO"		// "L" stands for linear interpolartion, "O" ordinate
	ydestwavname = ValidateName(ydest = 1)
	if(cmpstr(ydestWavname, ksValidateNameCanceled) == 0)
		return 0  // user canceled
	endif

	duplicate/O wav $ydestwavname/wave=ydestwav
	redimension/N=(-1,dimsize(pixel_new,0),-1,-1) ydestwav

	make/n=(dimsize(pixel_new,0))/O ydest_temp
	make/O/n=(dimsize(pixel_old,0)) wav_temp

	variable i,j,l
	l=0
	do
		j=0
		do
			i=0
			do
				// do interpolation

				wav_temp = wav[i][p][j][l]

				interpolate2/T=(interpmode)/I=3/E=2/Y=ydest_temp/X=pixel_new pixel_old, wav_temp

				ydestwav[i][][j][l] = ydest_temp[q]

				i+=1
			while(i < dims[0])

			j+=1
		while(j < dims[2])

		l+=1
	while(l < dims[3])
end

Threadsafe Function TSsin(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = y0 + A*sin(f*x + phi)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = offset
	//CurveFitDialog/ w[1] = A
	//CurveFitDialog/ w[2] = f
	//CurveFitDialog/ w[3] = phi


	return w[0] + w[1] * sin(w[2] * x + w[3])
End

Function skewedsin(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = y0 + A*sin(f*x + phi)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ w[0] = offset
	//CurveFitDialog/ w[1] = A
	//CurveFitDialog/ w[2] = f
	//CurveFitDialog/ w[3] = phi
	//CurveFitDialog/ w[3] = distortion


	return w[0] + w[1] * sin(w[2] * x + w[3] + triangle({0, w[4],w[2], w[3]},x))
End


//////////////////////////////


function TiltCorrectionPrep1(refwav)
	wave refwav
	newDataFolder/O/S $"TiltEvaluation"
	// - preparation of tempwav
	duplicate/O refwav, whitelight_rot
	VolumeTranspose(whitelight_rot, {1, 0, 2, 3}, dest = self())

	duplicate/O whitelight_rot whitelight_stripe_temp
	wavecrop(whitelight_stripe_temp, "[][0][0]", dest = self())

	display whitelight_stripe_temp
	ModifyGraph grid(bottom)=1,minor(bottom)=1,nticks(bottom)=10
	setAxis left 0,wavemax(refwav)
	MFCFinit()
	MFCFywavSelecProc("", 1, "whitelight_stripe_temp")
	MFCFxwavSelecProc("", 1, "_calculated_")

	print "Select Ywav and Xwav manually!"
	print "Then, run prep2..."
	print "\t",ProcsTemplate("TiltCorrectionPrep2","")
end

function TiltCorrectionPrep2([fixbg])
//function TiltCorrectionPrep2(NumPeaks, initoffset_pixel, separation_pixel,initHeight, initwidth, [fixbg])
	variable fixbg

	string funcname = "TiltCorrectionPrep2"
	string caller = setCaller("[" + funcname + "]")
	fixbg = selectnumber(paramisdefault(fixbg), fixbg, 1)
	DFREF df = $MFCFksDF()
	NVAR numALl = df:NumAll
	NVAR BLfit_Flag = df:BLFit_Flag
	wave W_coefALL = df:W_coefAll
	Wave W_coefnonzeroAll = df:W_coefnonzeroAll
	wave W_coefBg = df:W_coefBg

	wave wav = whitelight_stripe_temp
	variable pBegin = 0
	variable pEnd = dimsize(wav, 0)
	variable noiseLevel = 1e+2
	variable smoothingFactor = 3
	variable maxPeaks = 1e+2
	Variable peaksFound= AutoFindPeaks(wav,pBegin,pEnd,noiseLevel,smoothingFactor,maxPeaks)

	if(peaksFound > 0)
		wave W_AutoPeakInfo
		sort2D(WaveCrop(W_AutoPeakInfo,"[][0]", dest = free()), W_AutoPeakInfo, "", dest = self())
		numALL = peaksFound
		if(!fixbg)
			MFCFbuttonguessProc("")
	//		BLFit_Flag = 0
	//	 	W_coefBg[] = 0
	 	endif

		MFCFsetVarnumAllproc("",peaksFound,"","")

		wave W_AutoPeakInfo_cr = WaveCrop(W_AutoPeakInfo,"[][1]", dest = free(), verbose=-inf)
		W_coefAll[0][] = 0
		W_coefAll[1][] = W_AutoPeakInfo[q][2]
		W_coefAll[2][] = W_AutoPeakInfo[q][0]
		W_coefAll[3][] = median(W_AutoPeakInfo_cr)

		W_coefnonzeroAll[1][] = 1
		MFCFcalctrialCurves()
		MFCFdispfitcompproc("")
	endif
	LogInfo({"peaksFound: ", num2str(peaksFound)},getCaller(),kVerbose)
	LogInfo({"Repeat Prep2 unitl best initvalues are found"},getCaller(),kVerbose)
	LogInfo({"Then Run..."},getCaller(),kVerbose)
	LogInfo({ProcsTemplate("TiltCorrectionPrep3","")},getCaller(),kVerbose)
	setCaller("")
end

function TiltCorrectionPrep3([fixbg])
	variable fixbg

	fixbg = selectnumber(paramisdefault(fixbg), fixbg, 1)
	TiltCorrectionDofit(fixbg = fixbg)
	MFCFdispfitcompproc("")

	DFREF df = $MFCFksDF()
	duplicate/O df:W_coefAll W_coefAll_init

	// print "Inspect the fitting result"
	// print "If it looks good, then proceed to..."
	// print ProcsTemplate("TiltCorrectionPrep4","")
	// print "if multiple planes are to be fit,\r"
	// print "\t","concatenate {W_coefAll_init}, W_coefAll_initwaves"
	// print "then,"
	// print "\t","duplicate/O whitelight_rot whitelight_stripe_temp;wavecrop(whitelight_stripe_temp, \"[][0][n]\", dest = self()) // n = 1, 2, ..."
	// print "\t","Run Prep2 and Prep3 and\r"
	// Print "\t","concatenate {W_coefAll_init}, W_coefAll_initwaves"
	// print "After setting all the init paremters"
	// print ProcsTemplate("TiltCorrectionPrep4","") // with
end

function TiltCorrectionPrep4([suffix,fixbg, initwaves, extraTConstraints])
	string suffix
	variable fixbg
	wave initwaves
	string extraTConstraints

	variable timer = ticks
	fixbg = selectnumber(paramisdefault(fixbg), fixbg, 1)
	suffix = selectstring(paramisdefault(suffix), suffix, "")
	extraTConstraints = selectstring(paramisdefault(extraTConstraints), extraTConstraints, "")

	DFREF df = $MFCFksDF()
	wave whitelight_rot
	wave whitelight_stripe_temp
	wave W_coefAll_init
	wave W_coefALL = df:W_coefAll
	wave W_coefBg = df:W_coefBg
	make/n=4/FREE/D dims = dimsize(whitelight_rot, p)
	make/O/N=0 $("W_coefALLsav" + suffix)/wave=W_coefALLsav
	make/O/N=0 $("W_sigmaAllsav" + suffix)/Wave=W_sigmaAllsav
	wave W_sigma

	variable i, j
	j = 0
	do
		i = 0
		if(waveexists(initwaves))
			W_coefALL = initwaves[p][q][j]
		else
			whitelight_stripe_temp = whitelight_rot[p][i][j]
			TiltCorrectionPrep2()
			TiltCorrectionPrep3()
			wave W_coefAll_init
			W_coefALL = W_coefAll_init[p][q]
		endif
		make/O/N=0 W_coefALLsavtemp, W_sigmaAllsavtemp
		do
			whitelight_stripe_temp = whitelight_rot[p][i][j]
			wavestats/Q whitelight_stripe_temp
			if(V_npnts/dims[0] > .7)
				TiltCorrectionDofit(fixbg = fixbg, extraTConstraints = extraTConstraints)
				doUpdate
				concatenate {W_coefALL}, W_coefALLsavtemp
				deletepoints 0, 2, W_sigma
				concatenate {W_sigma}, W_sigmaAllsavtemp
			endif
			redimension/N=(dimsize(W_coefALL,0),dimsize(W_coefALL,1),i+1) W_coefALLsavtemp,W_sigmaAllsavtemp
			if(V_npnts/dims[0] < .7)
				W_coefALLsavtemp[][][i] = nan
			endif
			i += 1
		while(i < dims[1])
			if(j == 0)
				duplicate/O W_coefAllSavtemp, W_coefAllSav
				duplicate/O W_sigmaAllSavtemp, W_sigmaAllSav
			else
				concatenate/NP=1 {W_coefAllSavtemp}, W_coefAllSav
				concatenate/NP=1 {W_sigmaAllSavtemp}, W_sigmaAllSav
			endif
			LogInfo({"dimR: ", num2str(j)},"[TiltCorrectionPrep4]",1)
		j += 1
	while(j < dims[2])
	print "If the fitting finished without error, then proceed to..."
	print "\t",ProcsTemplate("TiltCorrectionPrep5","RS SERDSanalysis.ipf")
	print "Otherwise, additional constraint may help; list them as follows:"
	print "extraTConstraints = \"Kn < XX;\", K0,K1 for baseline, K2-K5 for 1st band,..."
	LogInfo({"Elapsed time: ", Secs2Time((ticks-timer)/60.15,5)}, "[TiltCorrectionPrep4]", kVerbose)
end

function TiltCorrectionPrep5([offset,polyorder])
	variable offset
	variable polyorder
	wave/Z mask
	polyorder = selectnumber(paramisdefault(polyorder), polyorder, 1)
	Loginfo({"polyorder: ", num2str(polyorder)}, "[TiltCorrectionPrep5]", 1)
	TiltCorrectionPrepPW(offset = offset)
	wave W_coefALLsav_peak,W_sigmaAllsav_peak
	if(WaveExists(mask))
		W_sigmaAllsav_peak *= mask[p]/mask[p]
	else
		duplicate/O W_sigmaAllsav_peak, mask
		mask = 1
	endif
	Approximate(W_coefALLsav_peak,W_sigmaAllsav_peak, polyorder)

	print "If the approximation looks ok, then proceed to..."
	print "\t",ProcsTemplate("TiltCorrectionPrep6","RS SERDSanalysis.ipf")
	print "In order to exclude selected points from the fitting, use 'mask' wav (0 or NaN elminates the point from fitting)"
	print "For closing all the fit windows... execute"
	print "\t","dowindows(winList(\"autofit_peak_comp*\",\";\",\"win:1\"),\"/K\")"
end

function TiltCorrectionPrepPW([offset])
	variable offset
	wave W_coefALLsav, W_sigmaAllsav

	wavecrop(W_coefALLsav, "[2][][]", dest = name("@_peak"))
	wavecrop(W_coefALLsav, "[3][][]", dest = name("@_width"))
	wave W_coefALLsav_peak, W_coefALLsav_width
	W_coefALLsav_peak += offset

	wavecrop(W_sigmaAllsav, "[2][][]", dest = name("@_peak"))
	wavecrop(W_sigmaAllsav, "[3][][]", dest = name("@_width"))
	wave W_sigmaAllsav_peak, W_sigmaAllsav_width

	sort2D(WaveCrop(W_coefALLsav_peak,"[][0]", dest = free()),W_sigmaAllsav_width,"", dest = self())
	sort2D(WaveCrop(W_coefALLsav_peak,"[][0]", dest = free()),W_coefALLsav_width,"", dest = self())
	sort2D(WaveCrop(W_coefALLsav_peak,"[][0]", dest = free()),W_sigmaAllsav_peak,"", dest = self())
	sort2D(WaveCrop(W_coefALLsav_peak,"[][0]", dest = free()),W_coefALLsav_peak,"", dest = self())
	matrixtranspose W_coefALLsav_peak
	matrixtranspose W_sigmaAllsav_peak
	matrixtranspose W_coefALLsav_width
	matrixtranspose W_sigmaAllsav_width
end

function TiltCorrectionPrep6([polyorder,discard])
	variable polyorder
	variable discard
	polyorder = selectnumber(paramisdefault(polyorder), polyorder, 1)
	discard = selectnumber(paramisdefault(discard), discard, 0)

	variable timer = ticks
	string funcname = "TiltCorrectionPrep6"
	string caller = "["+ funcname +"]"
	wave fit_W_coefALLsav_peak
	variable dim = dimsize(fit_W_coefALLsav_peak, 1)
	duplicate/O fit_W_coefALLsav_peak fit_W_coefALLsav_peak_cr
	fit_W_coefALLsav_peak_cr = nan
	fit_W_coefALLsav_peak_cr[][discard, dim - 1 - discard] = fit_W_coefALLsav_peak[p][q]
	wave fit_W_coefALLsav_peak = fit_W_coefALLsav_peak_cr
	duplicate/O fit_W_coefALLsav_peak, refpeaks
	wavecrop(refpeaks, "["+ num2str(floor(dimsize(fit_W_coefALLsav_peak,0)/2))+"][]", dest = self(), verbose = 0)
	LogInfo({"Reference Horizontal pixel:", "\t", num2str(floor(dimsize(fit_W_coefALLsav_peak,0)/2))}, "[TiltCorrectionPrep6]", 1)

	if(0)
		if(polyorder == 1)
			CurveFit/Q/ODR=2/N=1 line refpeaks[discard, dimsize(refpeaks,0) - 1 - discard]
		else
			CurveFit/Q/ODR=2/N=1 poly (polyorder+1), refpeaks[discard, dimsize(refpeaks,0) - 1 - discard]
		endif
		wave w_coef
		refpeaks = NaN
		refpeaks[discard, dim - 1 - discard] = poly(w_coef,p)
	else
		// refpeaks = NaN
		// refpeaks[discard, dim - 1 - discard] = WaveCrop(fit_W_coefALLsav_peak, "["+ num2str(floor(dimsize(fit_W_coefALLsav_peak,0)/2))+"][]", dest = free(), verbose = 0)[p]
	endif

	calcTiltCorrectionWlcaligparam(refpeaks, fit_W_coefALLsav_peak, polyorder, suffix = "_tilt")
	SVAR S_wavename
	wave destWav = $S_wavename
	note destWav, note(fit_W_coefALLsav_peak)
	setNoteProc(destWav, "TiltCorrectionPrep6", "polyorder", num2str(polyorder), newline = 1)
	setNoteProc(destWav, "TiltCorrectionPrep6", "discard", num2str(discard))
	LogInfo({"reference vertical pixel: ", varwave2list(refpeaks, listsepstr = ",")}, "[TiltCorrectionPrep6]", 1)
	LogInfo({"reference ragne (pixel): ", num2str(refpeaks[dim-1-discard] - refpeaks[discard]), "\t", "Separation (pixel): ", num2str(refpeaks[discard+1] - refpeaks[discard])}, "[TiltCorrectionPrep6]", 1)
end

static function Approximate(refwav,sigmawav, polyorder, [q1,q2, suffix])
	wave refwav
	wave sigmawav
	variable polyorder
	string suffix
	variable q1,q2

	variable verbose
	string funcname = "Approximate"
	string caller = setCaller("[" + funcname + "]")

	variable numPeak = dimsize(refwav,1)
	suffix = selectstring(paramisdefault(suffix),suffix,"")
	string destwavname = ValidateName("fit_" + nameofwave(refwav) + suffix,overwrite = 1)
	if(cmpstr(destWavname, ksValidateNameCanceled) == 0)
		setCaller("")
		return 0  // user canceled
	endif
	q1 = selectnumber(ParamIsDefault(q1),q1,0)
	q2 = selectnumber(ParamIsDefault(q2),q2,dimsize(refwav,0) - 1)
	variable dimx = dimsize(refwav,0)
	string win
	variable j

	dowindows(winList("autofit_peak_comp*",";","win:1"),"/K")
	duplicate/O refwav, $destwavname/wave=fit_refwav

	if(polyorder)
		for(j=0;j<numpeak;j+=1)
			if(polyorder == 1)
				CurveFit/Q/ODR=2/N=(!!verbose) line, refwav[q1,q2][j] /I=1/W=sigmawav[q1,q2][j]
//				CurveFit/Q/ODR=2/N=1 line, refwav[q1,q2][j]
				wave w_coef
				fit_refwav[][j] = w_coef[0] + w_coef[1] * p
			elseif(polyorder > 1)
				CurveFit/Q/W=0/N=(!!verbose) poly (polyorder+1),refwav[q1,q2][j] /I=1/W=sigmawav[q1,q2][j]
//				CurveFit/Q/W=0/N=1 poly (polyorder+1),refwav[q1,q2][j]
				wave w_coef
				fit_refwav[][j] = poly(W_coef,p)
			endif
			setNoteProc(fit_refwav,"Approxmate","w_coef",varwave2list(w_coef),newline=1)
		endfor
		setNoteProc(fit_refwav, "Approxmate", "polyorder", num2str(polyorder),newline = 1)

		matrixop/O res_M_temp_peakposi =  refwav - fit_refwav
		for(j=0;j<numpeak;j+=1)
			win = "autofit_peak_comp"+num2str(j)
			dowindow/F $win
			if(V_flag)
				dowindow/K $win
			endif
			display/K=1/N=$win refwav[][j],fit_refwav[][j]
			ErrorBars $nameofwave(refwav) SHADE= {0,0,(0,0,0,0),(0,0,0,0)},wave=(sigmawav[*][j],sigmawav[*][j])
			appendtograph/W=$win/L=Res_Left res_M_temp_peakposi[][j]
			setAxis left, wavemin(refwav,j*dimx,(j+1)*dimx-1),wavemax(refwav,j*dimx,(j+1)*dimx-1)
//			print wavemin(refwav,j*dims[1],(j+1)*dims[1]-1),wavemax(refwav,j*dims[1],(j+1)*dims[1]-1)
			ModifyGraph axisEnab(Left)={0,0.75}
			ModifyGraph axisEnab(Left)={0,0.75}
			ModifyGraph axisEnab(Res_Left)={0.8,1},freePos(Res_Left)=0,zero(Res_Left)=1
			ModifyGraph rgb($nameofwave(fit_refwav))=(0,0,65535)
			setaxis bottom, q1,q2
			setAxis/A=2 Res_left
			doupdate
		endfor
		tilerecentwindows(numpeak)
	endif
	string/g S_wavename = destwavname
	setCaller("")
end

static function calcTiltCorrectionWlcaligparam(peak_ref, fit_M_temp_peakposi, polyorder, [suffix])
	wave peak_ref
	wave fit_M_temp_peakposi
	variable polyorder
	string suffix

	suffix = selectstring(paramisdefault(suffix),suffix,"")
	string paramwavname = ValidateName("M_WlCalibParam" + suffix, overwrite = 1)
	if(cmpstr(paramwavname, ksValidateNameCanceled) == 0)
		return 0  // user canceled
	endif
	string paramsigmawavname = ValidateName("M_WlCalibParam_sigma" + suffix, overwrite = 1)
	if(cmpstr(paramsigmawavname, ksValidateNameCanceled) == 0)
		return 0  // user canceled
	endif

	variable dimx = dimsize(fit_M_temp_peakposi,0)
	make/n=(dimx,polyorder+1)/O $paramwavname/wave=M_WlCalibParam
	make/n=(dimx,polyorder+1)/O $paramsigmawavname/wave=M_WlCalibParam_sigma
	make/n=0/O W_chisq
	variable i
	variable V_fitOptions = 4
	variable V_FitError = 0
	duplicate/O peak_ref tempwav

	dowindow/F pause
	if(!V_flag)
		display/K=1/N=pause peak_ref vs tempwav
		ModifyGraph mode=3,marker=41
		tilerecentwindows(1)
	endif
	for(i=0;i<dimx;i+=1)
		tempwav = fit_M_temp_peakposi[i][p]
		if(polyorder == 1)
			CurveFit/Q/ODR=2/N=1 line peak_ref /X=tempwav/D/R
		else
			CurveFit/Q/ODR=2/N=1 poly (polyorder+1), peak_ref /X=tempwav/D/R
		endif
		ModifyGraph zero(Res_Left)=1
		wave w_coef,w_sigma
		wave reswav = $("res_" + nameofwave(peak_ref))
		reswav = selectnumber(numtype(tempwav[p]) == 2,reswav[p] ,nan)
		addtotail(W_chisq,V_chisq)
		V_FitError = 0

		M_WlCalibParam[i][] = w_coef[q]
		M_WlCalibParam_sigma[i][] = w_sigma[q]

		doUpdate
		if(0)
			pauseforuser pause
			display/K=1/N=pause peak_ref vs tempwav
			ModifyGraph mode=3,marker=41
			tilerecentwindows(1)
		endif
	endfor	// i
	string/g S_wavename = paramwavname
end

static function TiltCorrectionDofit([fixbg, extraTConstraints])
	variable fixbg
	string extraTConstraints
	extraTConstraints = selectstring(paramisdefault(extraTConstraints), extraTConstraints, "")
	DFREF df = $MFCFksDF()
	wave W_coefALL = df:W_coefAll
	wave W_coefholdAll = df:W_coefholdAll
	if(!fixbg)
		MFCFbuttonguessProc("")
	endif

	string/g root:Packages:RSMFCF:extraTConstraints = ""
	W_coefholdAll[2,3][] = 1
	W_coefholdAll[0][] = 1
	MFCFdofitProc("")

	W_coefholdAll[][] = p == 3
	W_coefholdAll[0][] = 1
	MFCFdofitProc("")

	string/g root:Packages:RSMFCF:extraTConstraints = extraTConstraints
	W_coefholdAll[][] = 0
	W_coefholdAll[0][] = 1
	MFCFdofitProc("")
end

function MakeWidthMapFromTiltEvaluation(srcWav)
	wave srcWav // template for retracting the image size
	TiltCorrectionPrepPW()
	wave W_coefALLsav_peak
	wave W_coefALLsav_width
	wave W_sigmaAllsav_width
	Wave xwav = summarize(W_coefALLsav_peak, "median", directionWav = {0}, dest = name("M_position_med"))

	MakeWidthMAP_worker(W_coefALLsav_width, xwav, {dimsize(srcWav, 0), dimsize(srcWav, 1)}, 2, sigmaWav = W_sigmaAllsav_width)
end

function makeWidthMAP_worker(srcwav, xwav, destWavSize, polyorder, [sigmaWav])
	wave srcWav
	wave xwav
	wave destWavSize
	variable polyorder
	wave/Z sigmaWav  // sigma for width

	make/n=4/FREE/D dims = dimsize(srcWav,p)
	variable i
	make/FREE/N=(destWavsize[0],destWavsize[1]) destWav

	if(!WaveExists(sigmaWav))
		duplicate/FREE srcwav, sigmaWav
		sigmawav = 1
	endif
	for(i=0; i < dims[0] ; i += 1)
		wave temp = wavecrop(srcWav, "[" + num2str(i) + "][]", dest = free(),verbose=0)
		wave sWav = wavecrop(sigmaWav, "[" + num2str(i) + "][]", dest = free(),verbose=0)
		if(polyorder == 1)
			curveFit/Q line, temp/X=xwav/I=1/W=sWav
			wave W_coef
			destWav[i][] = w_coef[0] + w_coef[1] *q
		elseif(polyorder > 1)
			curveFit/Q poly (polyorder+1), temp/X=xwav/I=1/W=sWav
			wave W_coef
			destWav[i][] = poly(w_coef, q)
		endif
	endfor

	duplicate/O destWav M_widthMAP
	//M_widthMAP[531,533][] = mean(wavecrop(M_widthMAP,"[" + num2str(p-2) +"," +num2str(p +2)+"][" + num2str(q) +"]", dest = free(),verbose=0))
end

function/WAVE getWidthConvKernel(widthMap,targetWidth)
	// worker function of ResolutionCorrection
	wave widthMap
	variable targetWidth

	make/n=4/FREE/D dims = dimsize(widthMap,p)
	make/FREE/N=(dims[0],dims[0]) M_kernel
	variable convWidth
	variable i

	for(i=0; i < dims[0] ; i += 1)
		convWidth = limit(sqrt(targetWidth^2 - widthMap[i]^2),1E-10,inf)
		M_kernel[][i] = gaussFWHMnorm({0,i,convWidth},p)
	endfor
	normalizeMD(M_kernel,{0}, "sum", dest = self())
	return M_kernel
end

function/wave ResolutionCorrection(srcWav, widthMap, [targetWidth, dest, verbose])
	// Apply resolution correction on images based on `widthMap` wave generated from `MakeWidthMap_workerFromTiltEvaluation`
	wave srcWav
	wave widthMap
	variable targetWidth
	wave/Z dest
	variable verbose

	variable timer = ticks
	string funcname = "ResolutionCorrection"
	string caller = setCaller("[" + funcname + "]")
	verbose = selectnumber(paramisdefault(verbose), verbose, kVerbose)

	targetWidth = selectNumber(ParamIsdefault(targetWidth), targetWidth, Wavemax(widthMap))
	make/n=4/FREE/D dims = dimsize(srcWav, p)
	variable i

	duplicate/FREE srcWav, destWav

	for(i=0; i < dims[0] ; i += 1)
		wave temp_src = wavecrop(srcWav, "[" + num2str(i) + "][]", dest = free())
		wave temp_map = wavecrop(widthMap, "[" + num2str(i) + "][]", dest = free())
		Wave M_kernel = getWidthConvKernel(temp_map,targetWidth)
		matrixOP/FREE temp_dest = M_kernel x temp_src
		destWav[i][][][] = temp_dest[q][r][s]
		LogProgress(i, dims[0], caller, verbose, step = floor(dims[0]/10))
	endfor

	note/k destwav, note(srcWav)
		setNoteHistory(destwav, funcname, "srcWav" ,getWavesDataFolder(srcWav,2), newline = 1)
		setNoteHistory(destwav, funcname, "widthMap" ,getWavesDataFolder(widthMap,2))
		if(!ParamIsdefault(targetWidth))
			setNoteHistory(destwav, funcname, "targetWidth" ,num2str(targetWidth))
		endif
		if(!ParamIsDefault(destWavname))
			setNoteHistory(destwav, funcname, "destWavname" ,destWavname, isStr = 1)
		endif
		if(!ParamIsDefault(overwrite))
			setNoteHistory(destwav, funcname, "overwrite" ,num2str(overwrite))
		endif

<<<<<<< Updated upstream:Igor User Procedures/RS slitconfocal_wlcalib_obsolete.ipf
	wave retWav = StdOutput2(destWav, PrepareDest(dest, srcWav, "", "_RC"))
=======
<<<<<<< Updated upstream:Igor User Procedures/RS slitconfocal_wlcalib.ipf
	wave retWav = StdOutput(srcWav, destWav, destWavname, overwrite)
=======
	wave retWav = StdOutput(destWav, PrepareDest(dest, srcWav, "", "_RC"))
>>>>>>> Stashed changes:Igor User Procedures/RS slitconfocal_wlcalib_obsolete.ipf
>>>>>>> Stashed changes:Igor User Procedures/RS slitconfocal_wlcalib.ipf
	LogInfo({"Elapsed time: ", Secs2Time((ticks-timer)/60.15,5)}, caller, verbose)
	setCaller("")
	return retWav
end
