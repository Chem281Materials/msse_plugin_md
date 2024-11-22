# Creating a Plugin Framework

As we've discussed, plugins can be a very powerful tool for enforcing code modularity, and for allowing external developers to easily modify or extend a code.
There are many different ways to create a framework for handling plugins, with every approach having its own set of benefits and drawbacks.
Ironically, highly sophisticated approaches aren't always ideal; it is common for developers to create an elaborate object-oriented plugin framework, only to be defeated by the [fragile base class problem](https://en.wikipedia.org/wiki/Fragile_base_class) or similar issues.

Sometimes, it is best for a plugin interface to be as minimalistic and flexible as possible when interacting with the core code, and this generally requires hard-coding as little specificity as possible in the interface itself.
For example, suppose we have a molecular dynamics code, and we want to support plugins that can perform the task of calculating the forces that are acting on the atoms.
This would make it easy to add support for new types of force fields, and would allow external users to easily create and distribute their own.
We might initially imagine that plugins should provide an `evaluate_forces` function with the following header:

```
extern "C"
void evaluate_forces(
       const int &nparticles,
       double &potential_energy,
       const double &box_size,
       const std::vector<std::array<double, 3>> &positions,
       std::vector<std::array<double, 3>> &forces)
```

The `evaluate_forces` function above has read access to the number of particles in the system (`nparticles`), the size of the simulation cell (`box_size`), and the nuclear coordinates of the atoms (`positions`).
It can also modify the forces (`forces`) and the potential energy of the system (`potential_energy`).
This might work fine for a while, but at some point someone might want to create a plugin where the `evaluate_forces` function requires additional information (for example, the charges on the atoms, or their element number).
This means that the `evaluate_forces` function for *all* plugins would need to be modified to accept additional arguments.
The above approach is highly inflexible - small modifications to the way in which the core code interacts with plugins will break existing plugins.

For a much more flexible alternative, consider a plugin function with the following header:

```
extern "C"
void evaluate_forces(
       std::map<std::string, std::shared_ptr<std::any>> &state)
```

The function accepts only a single argument, called `state`, which is a reference to a `map`.
The map's keys are strings, while its values are shared pointers that each point to an `std::any`.
This means that they could contain any type - some of the values in the map might ultimately correspond to integers, some to vectors, and others to custom classes.
In other words, our `state` argument could in principle contain *any and all* information we might need to communicate between the core code and the plugin.
For example, suppose `state` has a key called "nparticles", and this key corresponds to a shared pointer that manages an `std::any`, which in turn contains an `int`.
It would then be possible to extract a reference to this `int` by doing:

```
  int& nparticles = std::any_cast<int&>(*state["nparticles"]);
```

You could extract a `double`, `std::vector<std::array<double, 3>>`, or other type by following the same process.
If you at some point decide to pass additional information to the `evaluate_forces` function, you don't need to modify anything about existing plugins.
You can simply have the core code add the additional information to the `state` map that it sends to the plugins, and all plugins will have access to this information.


## Tasks

In this problem, you will develop a minimalistic plugin framework for a molecular dynamics code, following the strategy described above.
Much of the code for both the core executable and the plugin has already been written and is provided in this repository; your job is to implement all the elements required to establish the interface between the two.
The `executable` subdirectory contains both a `CMakeLists.txt` file and an `src/main.cpp` file.
The `src/main.cpp` file contains code that should be fairly familiar by now; it is a very simple molecular dynamics implementation.
There is one big change, however: all of the code associated with computing the forces between the atoms has been removed from `src/main.cpp`.

The `plugin` subdirectory contains an `src/plugin.cpp` file.
This file contains code for computing forces between atoms using the Lennard-Jones potential.
You'll be working to develop `src/plugin.cpp` into a proper plugin, and to modify the molecular dynamics executable to be able to use this plugin for evaluation of the forces.




### Task 1

Create a `CMakeLists.txt` file in the `plugin` directory that compiles the plugin.
The plugin does not need to be installable through CMake.

### Task 2

Modify the executable code (in `executable/src/main.cpp`) so that the end-user supplies the path to a plugin as a command-line argument.
The executable should then load the plugin (through `dlopen`).
Call the plugin's `initialize` function at the location in the code indicated by the comment `Call the plugin initialize function here`.
Call the plugin's `evaluate_forces` function at the location in the code with the comment `Call the plugin evaluate_forces function here`.

Properly calling the `initialize` and `evaluate_forces` functions will require that you prepare an `std::map` that includes all the data required by those function calls.
Creating this map will require modifying how you manage important data like the `forces` and `positions` vectors.
See the Hints section below for tips on how to do this cleanly and with minimal frustration.

### Task 3

In the plugin, write a template function having the following declaration:

```
template <typename T>
T& extract_from_state(std::map<std::string, std::shared_ptr<std::any>> &state, 
               const std::string key)
```

The `extract_from_state` function must return a reference of type `T` to the data managed by the `shared_ptr` in `state` that corresponds to the key "`key`".
In other words, this function is intended to provide a simple way to extract references to the values in `state`, without requiring you to strew `any_cast` all over your code. 
The `extract_from_state` function must check to confirm that `key` exists in `state`, and must throw an appropriate runtime error (with an informative error message) if it does not.
The `extract_from_state` function must also throw an appropriate runtime error (with an informative error message) in the event of a `bad_any_cast`.

### Task 4

In the plugin, modify the `evaluate_forces` forces function so that it calls `evaluate_lj_forces`.
If everything is set up correctly (in both your plugin and the executable), this should allow your executable to run simulations with a Lennard-Jones potential through the plugin.

### Task 5

Create a GitHub Action that runs whenever a push or PR is made to any branch.
The action should compile both the executable and the plugin, and run a simulation using the Lennard-Jones potential.

### Task 6

Write a summary (at least a couple paragraphs long) of the merits and limitations of this approach to handling plugins.
Try to anticipate possible problems that might arise over time for a project that follows this strategy.

### Task 7

A more common way to handle function calls with arbitrary numbers of arguments is to use [variadic function templates](https://www.geeksforgeeks.org/variadic-function-templates-c/).
Would it be possible to switch from our current strategy, in which we send a `std::map<std::string, std::shared_ptr<std::any>>` to our `evaluate_forces` function, to a strategy in which we use variadic function templates?
In other words, could we do something like the following?

```
extern "C"
template <typename... Types>
void evaluate_forces(Types... args) {
  ...
}
```

Explain why or why not.

## Hints

One major inconvenience that arises from storing important data in `std::shared_ptr<std::any>` types is that you might find yourself doing lots of awkward casts with `any_cast` whenever you want to access the data.
In the plugin, your `extract_from_state` function will help alleviate this problem.
In the context of your executable, a cleaner solution is to modify the `MDSimulation` class to maintain references to the data that is managed by each `shared_ptr` that you create.
If you do this correctly, you won't need to change much about the `main.cpp` code, aside from the class declaration and the `MDSimulation` initializer list.
For example, consider the following code, which manages the data for the number of particles within a `shared_ptr` named `nparticles_ptr`, while also maintaining a reference to the information that is managed by this smart pointer, which is named `nparticles`.
Throughout the rest of the code, `nparticles` can be used as though nothing has changed about the way this information is being managed.

```
class MDSimulation {
    ...
  private:
    std::shared_ptr<std::any> nparticles_ptr;
    int& nparticles; 
    ...
};
```

Note that references that are class members **must** be initialized in the class initializer list; you can't just declare a reference and then assign it later.
You can apply the same strategy to your `positions` and `forces` vectors, although you'll need to think carefully about how to handle the vectors within your class initializer list.
If *don't* do something along these lines, you'll likely end up needing to make lots of invasive changes to `main.cpp`.

## Answers
