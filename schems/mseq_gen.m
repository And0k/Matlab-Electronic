function [code]=mseq_gen(p, m, format) 
 
% Usage: [code] = mseq_gen(p, m, format) 
% 
% Parameters: 
%    p       = base number 
%    m       = number of shift registers 
%    format  = output numeric representation desired 
%                -- set to 0 for integer format (e.g. 0, 1 or 2) 
%                -- set to 1 for bipolar format (e.g. -1, 0 and 1)  
%                -- set to 2 for polyphase format (i.e. complex values) 
% Returns: 
%    code = binary (p=2) or non-binary (p>2) max-length shift register sequence, or 
%    code = [] if error occurs. 
% 
% Currently, p and m can take the following values: 
%       p          m 
%       2          2 ... 16 
%       3          2 ... 7 
%       5          2 ... 5 
%       7          2, 3, 4 
%       11         2, 3 
%       13         2, 3 
%       17         2, 3 
%       19         2, 3 
%       23         2, 3 
% 
% Written by Group 792, Aalborg University, Oct 2002 
 
 
% This mTable defines the combinations of p & m which are valid for the function: 
%    mTable[p, m] = 1,  if p and m is valid  
%                   0,  if p or m is not valid 
% 
%         1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 
mTable = [0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0;    % p = 1 
          0 1 1 1 1 1 1 1 1 1  1  1  1  1  1  1  0;    % p = 2 
          0 1 1 1 1 1 1 0 0 0  0  0  0  0  0  0  0;    % p = 3 
          0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0;    % p = 4 
          0 1 1 1 1 0 0 0 0 0  0  0  0  0  0  0  0;    % p = 5 
          0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0;    % p = 6 
          0 1 1 1 0 0 0 0 0 0  0  0  0  0  0  0  0;    % p = 7 
          0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0;    % p = 8 
          0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0;    % p = 9 
          0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0;    % p = 10 
          0 1 1 0 0 0 0 0 0 0  0  0  0  0  0  0  0;    % p = 11 
          0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0;    % p = 12 
          0 1 1 0 0 0 0 0 0 0  0  0  0  0  0  0  0;    % p = 13 
          0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0;    % p = 14 
          0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0;    % p = 15 
          0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0;    % p = 16 
          0 1 1 0 0 0 0 0 0 0  0  0  0  0  0  0  0;    % p = 17 
          0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0;    % p = 18 
          0 1 1 0 0 0 0 0 0 0  0  0  0  0  0  0  0;    % p = 19  
          0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0;    % p = 20 
          0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0;    % p = 21 
          0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0;    % p = 22 
          0 1 1 0 0 0 0 0 0 0  0  0  0  0  0  0  0];   % p = 23 
 
% Define primitive polynomials for different p: 
% The primitive polynomials' coefficiences are from  
% "Sequence Design for Communications Applications", Pingzhi Fan & Michael Darnell 
% ************ p = 2 ********************** 
%         h0 h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 h11 h12 h13 h14 h15 h16 
Xtable2= [1  1  0  0  0  0  0  0  0  0  0   0   0   0   0   0   0;  % m=1 
          1  1  1  0  0  0  0  0  0  0  0   0   0   0   0   0   0;  % m=2 
          1  1  0  1  0  0  0  0  0  0  0   0   0   0   0   0   0;  % m=3 
          1  1  0  0  1  0  0  0  0  0  0   0   0   0   0   0   0;  % m=4 
          1  0  1  0  0  1  0  0  0  0  0   0   0   0   0   0   0;  % m=5 
          1  1  0  0  0  0  1  0  0  0  0   0   0   0   0   0   0;  % m=6 
          1  1  0  0  0  0  0  1  0  0  0   0   0   0   0   0   0;  % m=7 
          1  0  1  1  1  0  0  0  1  0  0   0   0   0   0   0   0;  % m=8 
          1  0  0  0  1  0  0  0  0  1  0   0   0   0   0   0   0;  % m=9 
          1  0  0  1  0  0  0  0  0  0  1   0   0   0   0   0   0;  % m=10 
          1  0  1  0  0  0  0  0  0  0  0   1   0   0   0   0   0;  % m=11 
          1  1  0  0  1  0  1  0  0  0  0   0   1   0   0   0   0;  % m=12 
          1  1  0  1  1  0  0  0  0  0  0   0   0   1   0   0   0;  % m=13 
          1  1  0  1  0  1  0  0  0  0  0   0   0   0   1   0   0;  % m=14 
          1  1  0  0  0  0  0  0  0  0  0   0   0   0   0   1   0;  % m=15 
          1  1  0  0  0  0  0  0  0  0  0   0   0   0   0   0   1]; % m=16 
% ************ p = 3 ********************** 
%         h0 h1 h2 h3 h4 h5 h6 h7 
Xtable3=[ 1  1  0  0  0  0  0  0;   % m = 1 
          2  1  1  0  0  0  0  0;   % m = 2 
          1  2  0  1  0  0  0  0;   % m = 3 
          2  1  0  0  1  0  0  0;   % m = 4 
          1  2  0  0  0  1  0  0;   % m = 5 
          2  1  0  0  0  0  1  0;   % m = 6 
          2  0  1  0  0  0  0  1];  % m = 7 
% ************ p = 5 ********************** 
%         h0 h1 h2 h3 h4 h5 h6 h7 
Xtable5=[ 2  1  0  0  0  0  0  0;   % m = 1 
          2  1  1  0  0  0  0  0;   % m = 2 
          2  3  0  1  0  0  0  0;   % m = 3 
          2  2  1  0  1  0  0  0;   % m = 4 
          2  4  0  0  0  1  0  0];  % m = 5 
% ************ p = 7 ********************** 
%         h0 h1 h2 h3 h4 h5 h6 h7 
Xtable7=[ 2  1  0  0  0  0  0  0;   % m = 1 
          3  1  1  0  0  0  0  0;   % m = 2 
          2  6  0  1  0  0  0  0;   % m = 3 
          5  3  1  0  1  0  0  0];  % m = 4 
% ************ p = 11 ********************** 
%         h0 h1 h2 h3 h4 h5 h6 h7 
Xtable11=[0  0  0  0  0  0  0  0;   % m = 1 
          2  4  1  0  0  0  0  0;   % m = 2 
          4  1  0  1  0  0  0  0];  % m = 3 
% ************ p = 13 ********************** 
%         h0 h1 h2 h3 h4 h5 h6 h7 
Xtable13=[0  0  0  0  0  0  0  0;   % m = 1 
          2  1  1  0  0  0  0  0;   % m = 2 
          6  1  0  1  0  0  0  0];  % m = 3 
% ************ p = 17 ********************** 
%         h0 h1 h2 h3 h4 h5 h6 h7 
Xtable17=[0  0  0  0  0  0  0  0;   % m = 1 
          3  1  1  0  0  0  0  0;   % m = 2 
          3  1  0  1  0  0  0  0];  % m = 3 
% ************ p = 19 ********************** 
%         h0 h1 h2 h3 h4 h5 h6 h7 
Xtable19=[0  0  0  0  0  0  0  0;   % m = 1 
          2  1  1  0  0  0  0  0;   % m = 2 
          4  1  0  1  0  0  0  0];  % m = 3 
% ************ p = 23 ********************** 
%         h0 h1 h2 h3 h4 h5 h6 h7 
Xtable23=[0  0  0  0  0  0  0  0;   % m = 1 
          5  2  1  0  0  0  0  0;   % m = 2 
          3  1  0  1  0  0  0  0];  % m = 3 
 
% Check whether p or m is out of range: 
if (p > 23) || (p < 2) || (m > 17) || (m < 2) 
    % Error message and return empty code[] if error 
    error(['Value of p or m is not valid for the function! Type ' ...
'''help mseq_gen'' forfurther information!']); 
    code = [];
    return; 
end 
 
if (mTable(p, m) == 0) 
    % Error message and return empty code[] if error 
    error(['Value of p or m is not valid for the function! Type' ...
      '''help mseq_gen'' for further information!']); 
    code = []; 
    return; 
end 
 
% Calculate length of sequence 
N = p^m - 1; 
% Find the correct Xtable for certain p:
switch (p) 
case 2 
    Xtable = Xtable2; 
case 3 
    Xtable = Xtable3; 
case 5 
    Xtable = Xtable5; 
case 7 
    Xtable = Xtable7; 
case 11 
    Xtable = Xtable11; 
case 13 
    Xtable = Xtable13; 
case 17 
    Xtable = Xtable17; 
case 19 
    Xtable = Xtable19; 
case 23 
    Xtable = Xtable23; 
otherwise 
    % Error message and return empty code[] if error 
    error(['''Value of p is not valid for the function! Type ' ...
'''help mseq_gen'' for further information!']); 
    code = []; Xtable = [];
    return; 
end 
 
% Define the initial values of the sequence a0, a1, ... am-1 
IValue = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; 
% Imagine we have a Linear Feedback Shift Registers consists of m registers 
%    r1   r2   r3   ...  rm 
% Fill the initial state to registers 
register= IValue; 
 
% The first value of the sequence is the output of the registers 
%(or the value of the mth register) 
code= NaN(1,N);
code(1) = register(m); 
% Calculate the rest values of the sequence 
% based on initial values and the registers 
for i=2:N 
    % First, calculate the value of the first register value based on the feedback 
    %   temp = value of r1 = - h1*rm - h2*rm-1 - ... - hj*rm-j+1 - ... - hm*r1 
    temp = 0; 
    for j=1:m 
        temp = temp - Xtable(m, j)*register(m-j+1); 
    end 
    % Shifting the registers towards one step 
    register(2:m)= register(1:m-1); 
    % Assign the first register value 
    % Remember that the value is in GF(p) 
    register(1)  = mod(temp, p); 
    % Assign the output to code 
    code(i)= register(m); 
end 
 
% Convert the format to bipolar (if necessary) 
if format==1 
    code = round(code - (p-1)/2); 
end 
 
% Convert the format to polyphase (if necessary) 
if format==2 
    code=exp(sqrt(-1)*(2*pi*code)/p); 
end 