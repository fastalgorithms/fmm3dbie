
function selected = select_top_fraction(refinelist,errqlist,fraction)

if isempty(refinelist)
    selected = refinelist;
    return
end

% Support either one score per flagged patch or one score per mesh patch.
if numel(errqlist) == numel(refinelist)
    scores = errqlist(:);
elseif numel(errqlist) >= max(refinelist)
    scores = errqlist(refinelist);
    scores = scores(:);
else
    error(['errqlist is incompatible with refinelist: expected one ', ...
        'score per flagged patch or one score per mesh patch.'])
end

[~,order] = sort(scores,'descend');
nselect = max(1,ceil(fraction*numel(refinelist)));
selected = refinelist(order(1:nselect));

end