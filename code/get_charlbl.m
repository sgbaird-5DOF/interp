function charlbl = get_charlbl(n)
arguments
    n(1,1) double = 26
end
% GET_CHARLBL  get character labels, e.g. '(a)', '(b)', ... for figures.
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:n));
chars = chars.';
charlbl = strcat('(',chars,')');
end