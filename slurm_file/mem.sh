
while true; do
    ps -C mapmaking -o pid=,%mem=,vsz= >> /tmp/mem.log
    sleep 1
done 
