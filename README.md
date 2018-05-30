# Evolvable Hardware

This is the main repository for the C++ code needed to run the evolvable hardware experiments.

## What you'll need

Here's a list of things you'll need to run the program:

* [Quartus II Lite 16.1](http://dl.altera.com/16.1/?edition=lite) (version 17 and above seem to generate wrong circuits)
* A [DE-1 SoC Board](http://www.terasic.com.tw/cgi-bin/page/archive.pl?Language=English&No=836)
* A compiler able to generate ARM code. I used DS-5, which I recommend since the .cproject and .project files mean you can jump right into development.
* An SD Card with a Yocto Linux image. Further instructions can be found [here](https://rocketboards.org/foswiki/Documentation/AlteraSoCDevelopmentBoardYoctoGettingStarted#SoC_EDS)
* An Ethernet cable
* A terminal capable of making serial connections. On Windows, I used [TeraTerm](https://ttssh2.osdn.jp/index.html.en)
* A terminal capable of calling general Linux utilities (such as `ssh` and `scp`)

The other cables needed (miniusb and JTAG) should be available alongside your DE-1 board.

## What to do

* The first thing you need to do is compile the [circuit](https://github.com/VitorCBSB/circuito-genetico) that computes an individual's fitness score. Sit tight, this might take a while.
* After it's compiled, load it into your board using the Programmer tool inside Quartus II.
* Next, change the TesteLambda.cpp's parameters at the top of the file, after includes, to your desired circuit behavior. This is explained in more detail below.
* Compile it using your preferred method to an ARM binary.
* Make a serial connection to your DE1-SoC board.
* If your DE1-SoC board is connected to the internet (using the Ethernet cable), use `ifconfig` to find out what IP is assigned to it.
  * If it doesn't have an IP assigned to it even with the cable connected, call `udhcpc`.
* Send your compiled binary to it using `scp`.
* Run your binary.

## Parameters

This section explains what some of the parameters mean.

### On the board

The board switches can affect program execution in the following ways:

- SWs from 0 to 7 are a manual input to your program if SW 8 is up. If it's down, they do nothing.
- If SW 8 is down, SWs from 6 to 7 are used as 2 bit number to determine how many samples per sequence element.
  - 00 -> 100
  - 01 -> 500
  - 10 -> 1000
  - 11 -> 2000
- LEDRs from 0 to 6 indicate the output of the current individual being evaluated.
- LEDRs from 7 to 9 show the internal state machine's state.

### On the program

- `NUM_IN` is the desired circuit's number of inputs.
- `NUM_OUT` for number of outputs.
- `INITIAL_ROW_COUNT` is how many logic gates should the program start running with.
- `MUTATION_RATE` is a percentage (from 0 to 1) of how much of an individual should change during mutation.
- `LAMBDA` is the lambda parameter of the 1 + `lambda` method. This indicates how many mutated offsprings should be generated every generation.
- `MAX_GENERATIONS` is for how long should we look for a solution with the current amount of logic gates before increasing it by one.

An individual's behavior is specified by a sequence of triplets of (Input, Ouput, Validity) and can be seen in the `inputOutputValidSequences` function right below the above parameters.

Each element of the triplet is described as an 8 bit number. Output and Input are what a circuit should output when given this input. Validity is used to know which of the output bits of a circuit should be taken into account when computing its fitness. This is useful for circuits which don't have a valid output during the first few elements of the sequence (such as flip-flops).

