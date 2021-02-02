qcc -tags $file.c
sh page2html $file.c > $file.c.html
rm -r *.c.tags
# open -a chrome $file.c.html

