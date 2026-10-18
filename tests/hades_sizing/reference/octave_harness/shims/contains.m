function tf = contains(str, pat)
  if iscell(str)
    tf = cellfun(@(s) ~isempty(strfind(s, pat)), str);
  else
    tf = ~isempty(strfind(str, pat));
  end
end
