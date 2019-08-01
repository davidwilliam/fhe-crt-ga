# FHE
Ruby code for an FHE Library

## Requirements

This code requires Ruby installed on your system. There are [several options for downloading and installing Ruby](https://www.ruby-lang.org/en/downloads/ "Download Ruby").

This project uses only Ruby standard libraries, so once you have Ruby installed (version 2.2.3 and greater), you have everything required to run the code.

## Usage

### Running tests

Once Ruby is installed on your machine, from the command line and in the root directory of the project, run the tests to check the code with the following command:

`$ rake`

You should get a result similiar to the following:

```console
Run options: --seed 9109

# Running:

...........

Finished in 5.316182s, 2.0692 runs/s, 9.0290 assertions/s.

11 runs, 48 assertions, 0 failures, 0 errors, 0 skips
```

### Ruby Interactive Shell

You can also run code from the Ruby Interactive Shell (IRB). From the project's root directory, execute the following command on the terminal:

`$ irb`

You will see the IRB's prompt. Next, command snippets for specific cases that can be executed on IRB.

#### Key Generation

Require the file the will boot the entire project on IRB:

`> require './x'`

Create the 'x' object with the required secret and public variables by passing a configuration for depth (first argument) and security parameter (second argument):

`> x = X::FHE.new(8,128)`

Encryt the number 231:

`> c1 = x.encrypt(231)`

Encryt the number 209:

`> c2 = x.encrypt(209)`

Add c1 e c2:

`> c1_add_c2 = c1 + c2`

Multiply c1 e c2:

`> c1_mul_c2 = c1 * c2`

Decrypt c1_add_c2:

`> x.decrypt(c1_add_c2)`

As a result you should get:

`=> (440/1)`

Decrypt c1_mul_c2:

`> x.decrypt(c1_mul_c2)`

As a result you should get:

`=> (48279/1)`
