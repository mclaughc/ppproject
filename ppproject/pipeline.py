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
      print(stage)
      # Run the stage if it isn't complete yet.
      if (not stage.is_done()):
        stage.run()

      # Adjust done flag if needed.
      if (not stage.is_done()):
        is_done = False
      else:
        logger.debug("stage %s is done", stage.get_name())

    return is_done

  def run(self, num_iterations = 1):
    for i in range(0, num_iterations):
      # break out if we're finished early
      if (not self._run_iteration()):
        break

  def run_until_complete(self):
    done = False
    while (not done):
      done = self._run_iteration()
