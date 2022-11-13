Basically, install gitlab-runner and docker.
see here:
https://cylab.be/blog/30/running-gitlab-tests-locally-with-docker-ce
Then,

.. code-block:: bash

    gitlab-runner exec docker main


or, if using a locally build docker image

.. code-block:: bash

    gitlab-runner exec docker --docker-pull-policy=never main


See docs here https://docs.gitlab.com/runner/executors/docker.html#how-pull-policies-work
