function s = gm_sigma_from_ampvar(ampVar)
% GM_SIGMA_FROM_AMPVAR  Amplitude noise SD for gelmetricsMLE.m (R2 4a/4b, 6):
% robust SD (1.4826 * MAD) of the quiescent-segment amplitudes (analysis.m's
% ampVar), which estimates the spread of a fitted segment amplitude without a
% burst. Robust, because a few quiescent segments contain undetected bursts.
a = ampVar(isfinite(ampVar));
s = 1.4826 * median(abs(a - median(a)));
end
