rule install_cluster_3_0:
    output: "output/progs/bin/cluster-1.59/bin/cluster"
    threads: 1
    params: mem="10G"
    shell: """
    mkdir -p progs/bin
    cd progs/bin
    wget -N http://bonsai.hgc.jp/~mdehoon/software/cluster/cluster-1.59.tar.gz
    tar xvfz cluster-1.59.tar.gz
    cd cluster-1.59
    ./configure --without-x --prefix $PWD
    make 
    make install
"""