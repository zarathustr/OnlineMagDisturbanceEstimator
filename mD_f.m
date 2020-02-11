function mD = mD_f(mag)
global ax ay az mx my mz;

mD = ax * (mx + mag(1)) + ay * (my + mag(2)) + az * (mz + mag(3));
end