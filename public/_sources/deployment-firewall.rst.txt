Deployment using proxy
======================

In some networks it is necessary to configure a proxy for the primertool in order to access the internet. This can be
done by following the steps below:

1. Build the Docker image with the ENV variables for the proxy configuration. To do this, uncomment the lines below
inside the :code:`Dockerfile` and replace the values with your proxy information.

.. code-block:: bash

    ENV HTTP_PROXY="http://[USER]:[PASSWORD]@[PROXY_IP]:8080/"
    ENV HTTPS_PROXY="http://[USER]:[PASSWORD]@[PROXY_IP]:8080/"
    ENV NO_PROXY=localhost,[LOCAL_IP],0.0.0.0

This is important so that docker can access the internet to download the necessary dependencies of the image.

2. The files :code:`/etc/systemd/system/docker.service.d/http-proxy.conf` and :code:`~/.docker/config.json` (if they exist)
contain the proxy configuration for the docker service. If you are using a proxy, you must configure these files with
the proxy information. To do this, follow the steps below:

2.1. Create the file :code:`/etc/systemd/system/docker.service.d/http-proxy.conf` with the following content:

.. code-block:: bash

    [Service]
    Environment="HTTP_PROXY=http://[USER]:[PASSWORD]@[PROXY_IP]:8080/"
    Environment="HTTPS_PROXY=http://[USER]:[PASSWORD]@[PROXY_IP]:8080/"
    Environment="NO_PROXY=localhost,[LOCAL_IP],

2.2. Create the file :code:`~/.docker/config.json` with the following content:

.. code-block:: bash

    {
        "proxies":
        {
            "default":
            {
                "httpProxy": "http://[USER]:[PASSWORD]@[PROXY_IP]:8080/",
                "httpsProxy": "http://[USER]:[PASSWORD]@[PROXY_IP]:8080/",
                "noProxy": "localhost,[LOCAL_IP]"
            }
        }
    }

3. Restart the docker service:

.. code-block:: bash

    sudo systemctl daemon-reload
    sudo systemctl restart docker

4. Check if the proxy configuration is correct:

.. code-block:: bash

    sudo systemctl show --property=Environment docker

5. Run the Docker image with the following command:

.. code-block:: bash

    sudo docker run -it --net=host -p 8501:8501 primertool
