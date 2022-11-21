# Convert "e+" to numbers
# https://unix.stackexchange.com/questions/621256/convert-all-exponential-numbers-to-decimal-numbers-in-a-file-in-linux

awk '{ for (i=1;i<=NF;i++) if ($i+0 == $i && $i ~ /e/) $i = sprintf("%1.0f", $i) } 1' PCOS_UKBall.txt > PCOS_UKBall_m.txt
