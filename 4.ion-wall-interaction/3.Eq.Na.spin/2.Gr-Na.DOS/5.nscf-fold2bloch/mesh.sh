NGX=`grep 'dimension x,y,z NGX =' "../OUTCAR" | head -1 | awk '{split($0,array," ")} END{print array[5]}'`
NGY=`grep 'dimension x,y,z NGX =' "../OUTCAR" | head -1 | awk '{split($0,array," ")} END{print array[8]}'`
NGZ=`grep 'dimension x,y,z NGX =' "../OUTCAR" | head -1 | awk '{split($0,array," ")} END{print array[11]}'`
NGXF=`grep 'dimension x,y,z NGXF=' "../OUTCAR" | head -1 | awk '{split($0,array," ")} END{print array[4]}'`
NGYF=`grep 'dimension x,y,z NGXF=' "../OUTCAR" | head -1 | awk '{split($0,array," ")} END{print array[6]}'`
NGZF=`grep 'dimension x,y,z NGXF=' "../OUTCAR" | head -1 | awk '{split($0,array," ")} END{print array[8]}'`

echo "NGX=${NGX}"
echo "NGY=${NGY}"
echo "NGZ=${NGZ}"
echo "NGXF=${NGXF}"
echo "NGYF=${NGYF}"
echo "NGZF=${NGZF}"
