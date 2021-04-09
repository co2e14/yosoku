from ._anvil_designer import mainformTemplate
from anvil import *
import anvil.server
import anvil.media
import time

class mainform(mainformTemplate):
  def __init__(self, **properties):
    self.init_components(**properties)
    
  def asu_predict_click(self, **event_args):
    self.sg_in = self.spacegroup.selected_value
    self.uc_in = self.unitcell.text
    self.s_in = self.sequence.text
    self.s, self.asu_pred_tab, self.asu_pred_n_copies = anvil.server.call('matthewsrupp', self.sg_in, self.uc_in, self.s_in)
    self.predictiontable.text = self.asu_pred_tab
    self.asupermol.text = self.asu_pred_n_copies
    pass

  def phase_pred_click(self, **event_args):
    if self.s_in.isdigit():
      self.s = int(self.s_in)
    self.sg_in = self.spacegroup.selected_value
    self.uc_in = self.unitcell.text
    self.s_in = self.sequence.text
    self.d_min_in = self.highres.text
    self.d_min_in_static = self.highres.text
    self.asu_mol = int(self.asupermol.text)
    self.ref_per_s = anvil.server.call('predict', self.sg_in, self.uc_in, self.d_min_in, self.asu_mol, self.s)
    self.ref_per_s_static = self.ref_per_s
    self.fit_eq, self.res_vs_refl = anvil.server.call('resrange', self.sg_in, self.uc_in, self.asu_mol, self.s)
    self.predplot = anvil.server.call('makegraph', self.fit_eq, self.d_min_in_static, self.ref_per_s_static, self.res_vs_refl)
    self.plot_png.source = self.predplot
    pass

  def save_plot(self, **event_args):
    anvil.media.download(self.predplot)
    pass

  def reset(self, **event_args):
    self.spacegroup.selected_value = "P1"
    self.asupermol.text = 1
    self.unitcell.text = ""
    self.sequence.text = ""
    self.predictiontable.text = "Resetting..."
    time.sleep(0.6)
    self.predictiontable.text = "ASU prediction table will appear here..."
    self.highres.text = ""
    self.plot_png.source = "data:image/png;base64,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""
    pass
