from ResultsHandlerLib import *
import pylab as pl


res = ReadResults("results//p.txt")

times = [res.times[0], res.times[2], res.times[5], res.times[15], res.times[40], res.times[-1]]
for t in times:
	pl.plot(res.coord, res.getSolutionAtTime(t))

pl.grid(True)
pl.show()
