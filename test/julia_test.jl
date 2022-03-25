using Logging

io = open("log.txt", "w+")
logger = SimpleLogger(io)

with_logger(logger) do
    @info("a context specific log message")
end

@debug("a global msg")

close(io)