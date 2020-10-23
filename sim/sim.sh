g='N=1500, frq=0, up=1.0, psd=-1, times=1e3'; d="run/a03"; rm -rf $d/{cms,std/log}/*; mkdir -p $d
for i in {1000..1200}; do
    for m in {10..50..10}; do
        c="M=$m, $g, seed=$i"; s="saveRDS(r, \"{n:05d}.rds\")"
	echo "Rscript -e 'source(\"ini.R\"); r=sim($c, tu=2., naf=.00, key=\"TC0\", tag=\"N00\"); $s'"
	echo "Rscript -e 'source(\"ini.R\"); r=sim($c, tu=2., naf=.10, key=\"TC0\", tag=\"N10\"); $s'"
        echo "Rscript -e 'source(\"ini.R\"); r=sim($c, tu=1., naf=.00, key=\"TC1\", tag=\"N00\"); $s'"
	echo "Rscript -e 'source(\"ini.R\"); r=sim($c, tu=1., naf=.10, key=\"TC1\", tag=\"N10\"); $s'"
    done
done | hpcwp - -d$d -q20 -m8 -p4 --wtm 1 --ln '*.rds' --cp "*.R" --cp "R" --tag ${d##*/}

g='N=1500, naf=.1, frq=0, times=1e3'; d="run/a01"; rm -rf $d/{cms,std/log}/*; mkdir -p $d
for i in {1000..1200}; do
    for m in {10..30..5}; do
        c="M=$m, $g, seed=$i"; s="saveRDS(r, \"{n:05d}.rds\")"
        echo "Rscript -e 'source(\"ini.R\"); r=sim($c, up=1.0, psd=-1, L=NL, key=\"PD0\", tag=\"L=A\"); $s'"
	echo "Rscript -e 'source(\"ini.R\"); r=sim($c, up=1.0, psd=-1, L=$m, key=\"PD0\", tag=\"L=M\"); $s'"
	echo "Rscript -e 'source(\"ini.R\"); r=sim($c, up=.99, psd=NL, L=NL, key=\"PD1\", tag=\"L=A\"); $s'"
	echo "Rscript -e 'source(\"ini.R\"); r=sim($c, up=.99, psd=NL, L=$m, key=\"PD1\", tag=\"L=M\"); $s'"
    done
done | hpcwp - -d$d -q20 -m8 -p4 --wtm 2 --ln '*.rds' --cp "*.R" --cp "R" --tag ${d##*/}


g='N=1500, naf=0, frq=0, up=.99, psd=NL, times=1e3'; d="run/b01"; rm -rf $d/{cms,std/log}/*; mkdir -p $d
for i in {1000..1100}; do
    for m in {10..30..5}; do
        c="M=$m, $g, seed=$i"; s="saveRDS(r, \"{n:05d}.rds\")"
        echo "Rscript -e 'source(\"ini.R\"); r=sim($c, L=10, d=10, key=\"L10\", tag=\"D10\"); $s'"
	echo "Rscript -e 'source(\"ini.R\"); r=sim($c, L=20, d=10, key=\"L20\", tag=\"D10\"); $s'"
	echo "Rscript -e 'source(\"ini.R\"); r=sim($c, L=30, d=10, key=\"L30\", tag=\"D10\"); $s'"

	echo "Rscript -e 'source(\"ini.R\"); r=sim($c, L=10, d=20, key=\"L10\", tag=\"D20\"); $s'"
	echo "Rscript -e 'source(\"ini.R\"); r=sim($c, L=20, d=20, key=\"L20\", tag=\"D20\"); $s'"
	echo "Rscript -e 'source(\"ini.R\"); r=sim($c, L=30, d=20, key=\"L30\", tag=\"D20\"); $s'"

	echo "Rscript -e 'source(\"ini.R\"); r=sim($c, L=10, d=30, key=\"L10\", tag=\"D30\"); $s'"
	echo "Rscript -e 'source(\"ini.R\"); r=sim($c, L=20, d=30, key=\"L20\", tag=\"D30\"); $s'"
	echo "Rscript -e 'source(\"ini.R\"); r=sim($c, L=30, d=30, key=\"L30\", tag=\"D30\"); $s'"
    done
done | hpcwp - -d$d -q45 -m8 -p4 --wtm 2 --ln '*.rds' --cp "*.R" --cp "R" --tag ${d##*/}
