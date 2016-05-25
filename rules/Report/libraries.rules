rule download_librairies_files:
    output:
        "Report/lib/d3.min.js",
        "Report/lib/d3pie.min.js",
        "Report/lib/jquery-2.1.1.min.js",
        "Report/lib/jquery.cookie.min.js",
        "Report/lib/bootstrap.min.css",
        "Report/lib/bootstrap.min.js",
    shell:
        """
        wget https://webdav-r3lab.uni.lu/public/R3lab/IMP/static/d3.min.js -O Report/lib/d3.min.js
        wget https://webdav-r3lab.uni.lu/public/R3lab/IMP/static/d3pie.min.js -O Report/lib/d3pie.min.js
        wget https://webdav-r3lab.uni.lu/public/R3lab/IMP/static/jquery-2.1.1.min.js -O Report/lib/jquery-2.1.1.min.js
        wget https://webdav-r3lab.uni.lu/public/R3lab/IMP/static/js.cookie.js -O Report/lib/jquery.cookie.min.js
        wget https://webdav-r3lab.uni.lu/public/R3lab/IMP/static/bootstrap.min.css -O Report/lib/bootstrap.min.css
        wget https://webdav-r3lab.uni.lu/public/R3lab/IMP/static/bootstrap.min.js -O Report/lib/bootstrap.min.js
        """

rule local_librarie_files:
    output:
        'Report/lib/vizbin_logo.png',
        'Report/lib/IMP_logo.png',
        'Report/lib/imp.js',
        'Report/lib/imp.html'
    params:
        type = ''.join([t.lower() for t in TYPES])
    shell:
        """
        mkdir -p Report/lib
        cp {LIBDIR}/vizbin_logo.png {output[0]}
        cp {LIBDIR}/IMP_logo.png {output[1]}
        cp {LIBDIR}/imp.{params.type}.js {output[2]}
        cp {LIBDIR}/imp.{params.type}.html {output[3]}
        """
