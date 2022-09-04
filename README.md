profiner
======

Classy name? I tried my best. Hereby I present you: The implementation of our fine grained profiler for GPGPU kernels.

It allows you to record per-allocation cache rates by simulating and recording cache rates per memory reference and apply custom aggregation onto them. We evaluated our profiler against NVIDIA Nsight Compute. We refer to our paper for results and further information.

Build instructions
------
First, make sure that you have the [NVBit](https://github.com/NVlabs/NVBit) binaries (contents of the `nvbit/core` folder) in the `nvbit` folder. We used version 1.5.5 for our experiments. This would possibly look like this:
```
git clone [this repo] && cd [this repo]
cd /tmp
wget https://github.com/NVlabs/NVBit/releases/download/1.5.5/nvbit-Linux-x86_64-1.5.5.tar.bz2
tar -xvjf nvbit-Linux-x86_64-1.5.5.tar.bz2
cd -
cp -r /tmp/nvbit_release/core nvbit
```

Then, do it classically. Run CMake.
```
mkdir build && cd $_
cmake ..
make
```

How to run
------
```
PLUGINS=./profiler/plugin_bufferrates/libbufferrates.so LD_PRELOAD=./profiler/libmem_trace.so YOUR_CUDA_PROGRAM
```

License & Reference
------
Most parts of our program are licensed under the GPLv3 license included as LICENSE.md file.

If you decide to use our code or code based on this project in your application, please make sure to cite our pacific graphics paper:

```
@article{buelow2022finegrained,
    title = {Fine-Grained Memory Profiling of GPGPU Kernels},
    author = {von Buelow, Max and Guthe, Stefan and Fellner, Dieter W.},
    journal = {Computer Graphics Forum},
    issuetitle = {Pacific Graphics},
    volume = {41},
    number = {7},
    editor = {Umetani, Nobuyuki and Vouga, Etienne and Wojtan, Chris},
    date = {2022-10},
    publisher = {The Eurographics Association and John Wiley \& Sons Ltd.},
    venue = {Kyoto, Japan},
    pubstate = {forthcoming}
}
```

Contact
------
For any trouble with building, using or extending this software, please use the project's integrated issue tracker. We'll be happy to help you there or discuss feature requests.

For requests not matching the above, please contact the maintainer at maximilian.von_buelow(at)tu-darmstadt.de.
