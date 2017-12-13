try:
    from matplotlib.pyplot import plot, gca, ylim
    from matplotlib.collections import LineCollection, PatchCollection
    from matplotlib.patches import Rectangle
except:
    print('Failed to import matlab')
    def _st_plot(self,plot_lines=True,set_ylim=True,**kwargs):
        pass
else:
    def _st_plot(self,plot_lines=True,set_ylim=True,**kwargs):
        if plot_lines or set_ylim:
            N=max([sp.aff for sp in self])
        spikes = LineCollection([(\
        (self.timeBlockSize * sp.timeBlock + sp.time,sp.aff-0.5), \
        (self.timeBlockSize * sp.timeBlock + sp.time,sp.aff+0.5)) \
            for sp in self],**kwargs)
        gca().add_collection(spikes)
        if plot_lines:
            for i in range(N+2):
                gca().axhline(i+0.5)
        if set_ylim:
            ylim([-0.5, N+0.5])
