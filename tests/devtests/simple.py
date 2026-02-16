"""
Simple demonstration of STIsim.
"""
import stisim as sti

sim = sti.Sim(diseases=sti.HIV(), networks=sti.StructuredSexual())
sim.run()
sim.plot()