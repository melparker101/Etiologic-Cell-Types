# This command trims the rsID column of the premunge sumstats file
# i.e. rs140052487:C:A to rs140052487

sed 's/:[^ ]*//g' bmi.giant-ukbb.txt > bmi_check_premunge.txt
