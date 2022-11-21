# For dichotomous traits, change the column name of Neff to N, as CELLECT does not recognise Neff
sed '0,/Neff/{s/Neff/N/}' PCOS_premunge2.txt > PCOS_premunge_N.txt
