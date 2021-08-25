⍝ Get a free trial version of dyalog from https://www.dyalog.com/
⍝ Works with version 18.0.40684
∇R←SP DMN_TIES c;m;n;TIES;k;diag;u;v;diag_bit;d
(m n)←⍴¨SP
SP←SP[⍋SP←∊SP]
TIES←(¯1↓SP≠1⌽SP),1
⎕←'Ties=',TIES
k←1
diag←1
u←0
⎕←'m n =',m n
⍝ LOOP:
:Repeat
  ⎕←'k=',k
  u←u,1+¯1↑u
  v←k-u
  diag_bit←(u≤m)∧(v≤n)∧(u≥0)∧v≥0
  ⍝(u v)←diag_bit/¨u v
  u←diag_bit/u
  v←diag_bit/v
  ⎕←'u=',u
  ⎕←'v=',v
  d←|(u÷m)-v÷n
  diag←diag_bit/(diag,0)+0,diag
  ⎕←'diag=',diag
  diag←∊(1 0=TIES[k])/(diag×c>d)(diag)
  ⎕←'diag=',diag
  k←k+1
:Until  (m+n)<k
⍝ →LOOP×~(m+n)≥k←k+1
⎕←'final diag=',diag
R←1-diag÷m!m+n
∇
S←(1 2 3 4 5)(6 7 8 9 10 11 12)
prob←S DMN_TIES 3÷7
prob
S←(1 2 2 3 3)(1 2 3 3 4 5 6)
prob←S DMN_TIES 3÷7
prob
