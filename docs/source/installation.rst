.. _getting-started:


install
=======

The easiest way to install ``abutils`` locally (on macOS or Linux) is to use ``pip``::

    $ pip install abutils

``abutils`` does not run natively on Windows, but Windows users can use Docker_ 
(the brineylab datascience_ Docker container contains the entire ab[x] toolkit,
which includes ``abutils``)::

    $ docker pull brineylab/datascience
    $ docker run -it brineylab/datascience

Stable_ and development_ versions of ``abutils`` can also be downloaded from GitHub. 
You can manually install the latest version of ``abutils`` with::

    $ git clone https://github.com/briney/abutils
    $ cd abutils/
    $ python setup.py install


.. _Docker: https://www.docker.com/
.. _datascience: https://hub.docker.com/repository/docker/brineylab/datascience/general
.. _stable: https://github.com/brineylab/abutils/releases
.. _development: https://github.com/briney/abutils

