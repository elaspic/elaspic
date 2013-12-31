prefix=../input/configFile-unique-
suffix=.ini

for i in {2,8,9,12}
do
submitjob 72 -m 8 -c 2 python pipeline.py -c $prefix$i$suffix
done
