awk -v OFS='\t' '
    BEGIN {print "low","high","concordant","discordant","rate","s.e.m.","expected_rate"}
    {
        expected_rate=$10;
        sign=$7;
        for (i = 0; i < 1; i+=0.05) {
            bin = i OFS i+0.05;
            if (expected_rate > i && expected_rate <= i + 0.05) {
                if ( sign > 0) {
                    concordant[bin] += 1;
                } else {
                    discordant[bin] += 1;
                }
                expected[bin] += expected_rate;
            }
        }
    }
    END {
        for (bin in concordant) {
            c = concordant[bin]+0;
            d = discordant[bin]+0;
            x = expected[bin]/(c+d);
            n = c+d;
            p = c/n;
            print bin,c,d,p,sqrt(p*(1-p)/n),x;
        }
    }
' | sort -k1,1n | column -t
