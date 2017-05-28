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
	@echo 'Invoking: GCC C++ Compiler 4 [arm-linux-gnueabihf]'
	arm-linux-gnueabihf-g++ -std=c++1y -O1 -finline-functions -fstrict-aliasing -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


