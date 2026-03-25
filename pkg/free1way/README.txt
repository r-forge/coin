
making changes:

- change code in free1way.w
- uncomment .Rbuildignore
- make / make all
- build / check

to stats

- NAMESPACE
- remove stats::: for C_dpermdist2’
- cp *R
- \link[free1way] -> \link[stats] in ppplot.Rd
- cp *Rd
- make
- make check-all
- update stats-Ex.Rout.save
- svn diff > patch

