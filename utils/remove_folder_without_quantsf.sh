for i in ../hERV_Work/results/salmon/*; do
    if [ ! -e $i/quant.sf ]; then
        rm -r $i
    fi
done