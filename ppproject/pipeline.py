import logging

logger = logging.getLogger("Pipeline")

class Pipeline:
  def __init__(self):
    self.stages = []

  def add_stage(self, stage):
    if (not stage.is_ready()):
      stage.make_ready()

    logger.debug("added stage '%s'", stage.get_name())
    self.stages.append(stage)

  def _run_iteration(self):
    # Start with a done flag of true, we'll clear it if anything isn't done.
    is_done = True
    for stage in self.stages:
      # Run the stage, and adjust done flag if needed.
      stage.run()
      if (not stage.is_done()):
        is_done = False
      else:
        logger.debug("stage %s is done", stage.get_name())

    return is_done

  def run(self, num_iterations = 1):
    for i in range(0, num_iterations):
      # break out if we're finished early
      if (self._run_iteration()):
        break

  def run_until_complete(self):
    done = False
    while (not done):
      done = self._run_iteration()


class PipelineStage:
  def __init__(self, name = ""):
    """Creates a new pipeline stage, with the specified name. If the name is empty,
    a name based on the string representation of the object will be generated."""
    # Give a name based on the pointer if one wasn't specified.
    if (len(name) == 0):
      self.name = str.format("filter_%s" % (repr(self)))
    else:
      self.name = name

  def get_name(self):
    """Returns the name of this pipeline stage."""
    return self.name
    
  def make_ready(self):
    """Performs any late initialization required by the stage."""
    pass

  def is_ready(self):
    """Returns true if the stage is ready to start executing."""
    return True

  def is_done(self):
    """Returns true if at the current point in time, this stage is completed."""
    return True

  def run(self):
    """Execute this pipeline stage, processing any pending input."""
    pass
