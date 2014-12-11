rm -f printable.rst
echo "#############" >  printable.rst
echo "Documentation" >> printable.rst
echo "#############" >> printable.rst
echo "" >> printable.rst
cat cont_introduction.rst >> printable.rst
echo "" >> printable.rst
cat cont_codeintroduction.rst >> printable.rst
echo "" >> printable.rst
cat cont_testoperators.rst >> printable.rst
echo "" >> printable.rst
cat cont_nondim.rst >> printable.rst
echo "" >> printable.rst
cat cont_simple.rst >> printable.rst
echo "" >> printable.rst
cat cont_solution.rst >> printable.rst
echo "" >> printable.rst
echo "######################" >> printable.rst
echo "Appendix 1: Test codes" >> printable.rst
echo "######################" >> printable.rst
echo "" >> printable.rst
cat cont_smith_hutton_script.rst >> printable.rst
echo "" >> printable.rst
cat cont_smith_hutton_diff_script.rst >> printable.rst
echo "" >> printable.rst
cat cont_inclined_diff_script.rst >> printable.rst
echo "" >> printable.rst
echo "#########################" >> printable.rst
echo "Appendix 2: Solution code" >> printable.rst
echo "#########################" >> printable.rst
echo "" >> printable.rst
cat cont_square_cylinder_script.rst >> printable.rst
echo "" >> printable.rst
echo "#################################" >> printable.rst
echo "Appendix 3: Library documentation" >> printable.rst
echo "#################################" >> printable.rst
echo "" >> printable.rst
for f in pyCFD*.rst; do cat $f >> printable.rst; echo "" >> printable.rst; done
echo "" >> printable.rst
# cat modules.rst >> printable.rst
vim printable.rst -s make_printable.vim_script
