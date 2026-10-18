function c = splitlines(txt)
  c = strsplit(txt, {"\r\n", "\n", "\r"}, 'CollapseDelimiters', false);
  c = c(:);
end
