function colSum(arr) {
    return arr.reduce((a, b) => a.map((x, i) => x + b[i]))
}

function computeCCF(V, T, omega) {
    n_samples = V.length;
    var ccfs = Array(n_samples).fill(0.0);
    for(let i = 0; i < n_samples; i++)
    {
        if((V[i] > 0) && (T[i] > 0) && (omega[i] > 1e-3))
        {
            ccfs[i] = (V[i]/(T[i]*omega[i]));
        }
    }
    ccfs = minimum(ccfs, Array(n_samples).fill(1.0))
    return ccfs;
}

function computeCWACCFs(c, supervariants, F)
{
    let nSamples = F[0].length;
    let numer = Array(nSamples).fill(0);
    let denom = Array(nSamples).fill(0);
    for(let i = 0; i < c.length; i++)
    {
        let id = c[i];
        for(let j = 0; j < nSamples; j++)
        {
            let fIndex = Number(id.slice(1))
            numer[j] += F[fIndex][j]*supervariants[id][totalReads][j];
            denom[j] += supervariants[id][totalReads][j];
        }
    }
    var CWACCFs = [];
    for(let j = 0; j < nSamples; j++)
    {
        CWACCFs.push(Math.max(1e-30,Math.min(1-1e-3, numer[j]/denom[j])));
    }
    return CWACCFs;
}

function ccfError(c, supervariants, treeCCFs)
{
    var err = 0;
    let nSamples = treeCCFs.length;
    for(let id of c)
    {
        ccfs = Array(nSamples).fill(0);
        for(let i = 0; i < nSamples; i++)
        {
            let V = supervariants[id][varReads][i];
            let T = supervariants[id][totalReads][i];
            let omega = supervariants[id][varReadProb][i];
            if((V > 0) && (T > 0) && (omega > 0))
            {
                ccfs[i] = V/(T*omega);
            }
        }
        ccfs = minimum(ccfs, Array(n_samples).fill(1.0));
        err += math.abs(math.subtract(ccfs, treeCCFs)).reduce((partialSum, a) => partialSum + a, 0);
    }
    return err;
}

function submitFiles(e) {  
    e.preventDefault();
    let ssmfile = e.target.ssm.files[0];
    let paramsfile = e.target.params.files[0];
    let treefile = e.target.tree.files[0];
        
    let promises = [];

    if(paramsfile)
    {
        // promise for params file
        let filePromise1 = new Promise(resolve => {
            // load params file 
            var reader1 = new FileReader();
            reader1.onload = function() {
                resolve(JSON.parse(reader1.result));
            };
            reader1.onerror = function() {
                console.log(reader1.error);
                resolve({});
            };
            reader1.readAsText(paramsfile);
        });
        promises.push(filePromise1);
    }

    if(ssmfile)
    {
        // promise for ssm file
        let filePromise2 = new Promise(resolve => {
            // load ssm file
            var reader2 = new FileReader();
            reader2.onload = function() {
                resolve(ssmToVariants(reader2.result));
            };
            reader2.onerror = function() {
                console.log(reader2.error);
                resolve({});
            };
            reader2.readAsText(ssmfile);
        });
        promises.push(filePromise2);
    }

    if(treefile)
    {
        // promise for npz file
        let filePromise3 = new Promise(resolve => {
            var reader3 = new FileReader();
            reader3.onload = function() {
                resolve(JSON.parse(reader3.result));
            };
            reader3.onerror = function() {
                console.log(reader3.error);
                resolve({});
            };
            reader3.readAsText(treefile);
        });
        promises.push(filePromise3);

    }

    if(promises.length > 0)
    {
        Promise.all(promises).then(values => {
            const params = values[0];
            const [variants, n_samples] = values[1];
            const treeData = values[2];
            var supervariants = makeSVs(variants, params[clustersKey], n_samples);
            if(params)
            {
                localStorage.setObject(paramsKey, params);
            }
            if(variants)
            {
                localStorage.setObject(variantsKey, variants);
            }

            if(supervariants)
            {
                localStorage.setObject(supervariantsKey, supervariants);
            }
            
            if(treeData)
            {
                localStorage.setObject(parentsKey, treeData[parentsKey]);
            }

            if(treeData)
            {
                localStorage.setObject(FKey, treeData[FKey]);
            }

            localStorage.setObject(useNameKey, false);
            localStorage.setObject(treeIndexKey, 0);
            
            let event = new CustomEvent("dataLoaded");
              
            window.dispatchEvent(event);

        });
    }
  }

function ssmToVariants (tsv) {
    var data = {};
    var vid, n_samples;
    rows = tsv.split("\n").slice(1);

    for(let row of rows) {
        var rowData = row.split("\t");

        if (rowData[0])
        {
            vid = rowData[0];
            data[vid] = {
                name:rowData[1],
                var_reads:rowData[2].split(",").map(Number),
                total_reads:rowData[3].split(",").map(Number),
                var_read_prob:rowData[4].split(",").map(Number)
            };
        }
    }

    n_samples = data[Object.keys(data)[0]][varReads].length;
    return [data, n_samples];
};

function convertVarReads(V, T, omega)
{
    let V_hat = [...V];
    let N_hat = math.dotMultiply(T.map((x) => 2*x), omega).map((x) => Math.round(x));
    V_hat = minimum(N_hat, V_hat);
    return V_hat;
}

function convertTotalReads(T, omega)
{
    let N_hat = math.dotMultiply(T.map((x) => 2*x), omega).map((x) => Math.round(x));
    return N_hat;
}

function makeSVs(variants, clusters, n_samples)
{
    supervariants = {};
    for(let i = 0; i < clusters.length; i++)
    {
        c = clusters[i];
        let id = "S" + i;
        let nodeName = id;

        if(c.length == 1) {
            nodeName = variants[c[0]][nameKey];
            console.log(nodeName);
        } 

        supervariants[id] = {
            id:id,
            name:nodeName,
            vars:c,
            vars_names:c.map((vid,_) => variants[vid][nameKey]) ,
            var_reads:colSum(c.map((vid, _) => convertVarReads(variants[vid][varReads], variants[vid][totalReads], variants[vid][varReadProb]))),
            total_reads:colSum(c.map((vid, _) => convertTotalReads(variants[vid][totalReads], variants[vid][varReadProb]))),
            var_read_prob:Array(n_samples).fill(0.5)
        }
    }
    return supervariants;
}