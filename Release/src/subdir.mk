################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/TesteLambda.cpp 

OBJS += \
./src/TesteLambda.o 

CPP_DEPS += \
./src/TesteLambda.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: ARM C++ Compiler 6'
	armclang --target=aarch64-arm-none-eabi -xc++ -fno-exceptions -O2 -MD -MP -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


