function s = gm_scorer()
% GM_SCORER  Returns the scorer selected by the GM_SCORER environment
% variable ('orig' when unset). Used by gm_score.m and to tag check outputs.
s = getenv('GM_SCORER');
if isempty(s), s = 'orig'; end
end
