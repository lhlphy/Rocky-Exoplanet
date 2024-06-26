mkdir ~/.ssh -p
echo "ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABgQDew2o3f04gAH4xYmIadw7j7tP2dt8d9LgjNJ9VClAa0DCZ7tDdATUuYLaW6X8ZTaQ9Sv7uS/oM4GvGV0CzPiOOWyEddrAtB2rFuX5MVISHNwJHp/Vff5J7elMB4j4teklCsatsO0Q6HWiHaMCW68wdsJeoz2xBvdRFcd9pAScZWUVPXSboJpf5CiAYUewiSY+cmi1oaQ7jwVuIQ+Y3s8RYC/orQKDgAtniNR/EXrB6LpnZWi3WQOht++XdV9mrJ6ZGvr96EO+aIge3qKu/E5jxFy45Id2+9ThjzGW4i/li9DmOWBELRuQpgvYAYHvyJQD8rZLXh+tJY9GU3cetNxzMNIyFVnhGRgJH3hOwXj46JeKLpovI4fqa+qVxAddI5AsIE28HBsQ+EzFArXbxWZ9fw5M5Tq3CpkehLCNWMhrTh6ymOftxSJIGejwoDgnw7TRWZlKRJ8LSZuh4MQjbqvTD02cM2CAK/1c0xvI8f2CHLQDx6igoQd0Ls5lkCRLYKr8= dell@DESKTOP-KLE2D52" > ~/.ssh/authorized_keys
git config --global user.name "Lignting2"
git config --global user.email "lhlphy@qq.com"
# 设置代理，http 和 https 都加上代理，代理到 http://127.0.0.1:1087 这个 vpn 本地 IP
git config --global http.proxy http://127.0.0.1:1087 
git config --global https.proxy http://127.0.0.1:1087

# 取消代理
git config --global --unset http.proxy
git config --global --unset https.proxy

# 查看代理
git config --global --get http.proxy
git config --global --get https.proxy

# 查看全局所有配置
git config --global --list

# 设置代理
git config --local http.proxy http://127.0.0.1:1087 
git config --local https.proxy http://127.0.0.1:1087

# 取消代理
git config --local --unset http.proxy
git config --local --unset https.proxy