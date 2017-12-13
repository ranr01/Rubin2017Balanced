TARGET ?= _ST.so

SRCS := $(shell ls *.cpp)
OBJS := $(addsuffix .o,$(basename $(SRCS)))
DEPS := $(OBJS:.o=.d)

INC_DIRS := ../headers /usr/include/python3.5m
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS ?= -O3 -DCPP11 -DNDBUG -DUSE_PYTHON $(INC_FLAGS) -std=c++11 -fPIC -MMD -MP
LDFLAGS ?= -lboost_python-py35 -shared -fPIC

$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS)

.PHONY: clean
clean:
	$(RM) $(TARGET) $(OBJS) $(DEPS)

-include $(DEPS)



