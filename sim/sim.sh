g='N=1500, M=2, psd=-1, times=1e3'; d="run/a03"; rm -rf $d/{cms,std/log}/*; mkdir -p $d
for i in {1000..1200}; do
    for l in {5..20..5}; do
        c="L=$l, $g, seed=$i"; s="saveRDS(r, \"{n:05d}.rds\")"
	echo "Rscript -e 'source(\"ini.R\"); r=sim($c, tu=2., naf=.00, key=\"TC0\", tag=\"N00\"); $s'"
	echo "Rscript -e 'source(\"ini.R\"); r=sim($c, tu=2., naf=.10, key=\"TC0\", tag=\"N10\"); $s'"
        echo "Rscript -e 'source(\"ini.R\"); r=sim($c, tu=1., naf=.00, key=\"TC1\", tag=\"N00\"); $s'"
	echo "Rscript -e 'source(\"ini.R\"); r=sim($c, tu=1., naf=.10, key=\"TC1\", tag=\"N10\"); $s'"
    done
done | hpcwp - -d$d -q20 -m8 -p4 --wtm 1 --ln '*.rds' --cp "*.R" --cp "R" --tag ${d##*/}

g='N=1500, M=1, b=1, e=4, kgp.psd=-1, kgp.ucr=2, tol.cor=.05, tol.egv=1e-8, times=2e2'; d="run/a00"; rm -rf $d*; mkdir -p $d
for i in {1000..1200}; do
    for l in 10*{1..5}; do
        c="$g, L=$l, seed=$i"; s="saveRDS(r, \"{n:05d}.rds\")"
        echo "Rscript -e 'source(\"ini.R\"); r=sm2($c, a=.0, na=.0, key=\"NA0\", tag=\"a=0\"); $s'"
        echo "Rscript -e 'source(\"ini.R\"); r=sm2($c, a=.1, na=.0, key=\"NA0\", tag=\"a=1\"); $s'"
	echo "Rscript -e 'source(\"ini.R\"); r=sm2($c, a=.0, na=.1, key=\"NA1\", tag=\"a=0\"); $s'"
        echo "Rscript -e 'source(\"ini.R\"); r=sm2($c, a=.1, na=.1, key=\"NA1\", tag=\"a=1\"); $s'"
    done
done | hpcwp - -d$d -q20 -m8 -p4 --wtm 1 --ln '*.rds' --cp "*.R" --cp "R" --tag ${d##*/}

for f in run/a??; do $f/sub.sh; done

for f in run/a??; do Rscript -e "source('rpt.R'); plt.pow('$f')"; done
